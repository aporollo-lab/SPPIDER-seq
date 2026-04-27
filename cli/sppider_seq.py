import argparse
import os
import random
import re
import time
import urllib.request
import zipfile
from glob import glob
from io import StringIO
from itertools import product

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import torch
    import torch.nn as nn
    from Bio import SeqIO
    from scipy.stats import norm
    from statsmodels.stats.multitest import fdrcorrection
    from transformers import EsmModel, EsmTokenizer
    IMPORT_ERROR = None
except ModuleNotFoundError as exc:
    plt = None
    np = None
    torch = None
    nn = None
    SeqIO = None
    norm = None
    fdrcorrection = None
    EsmModel = None
    EsmTokenizer = None
    IMPORT_ERROR = exc


ESM_MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
PEPTIDE_MODEL_URL = "https://raw.githubusercontent.com/aporollo-lab/SPPIDER-seq/main/models/crossattn_pep_run07_best.pt"
RECEPTOR_MODEL_URL = "https://raw.githubusercontent.com/aporollo-lab/SPPIDER-seq/main/models/crossattn_rec_run27_best.pt"
PEPTIDE_MODEL_PATH = "models/peptide_model.pt"
RECEPTOR_MODEL_PATH = "models/receptor_model.pt"
PEPTIDE_MODEL_FALLBACK = "models/crossattn_pep_run07_best.pt"
RECEPTOR_MODEL_FALLBACK = "models/crossattn_rec_run27_best.pt"

NUM_HEADS = 16
MAX_TOKENS = 1024
STRIDE = 512
PRED_CUTOFF = 0.50
EMBEDDING_CACHE_MAX_SEQS = 256
EMBEDDING_CACHE_MAX_LEN = 3000

device = None
tokenizer = None
esm_model = None
EMBED_DIM = 1280

embedding_cache = {}
peptide_site_model = None
receptor_site_model = None

use_null_background = False
show_plots = False
pred_time_total_s = 0.0
pred_time_pairs = {}
pred_num_pairs = 0
base_scramble_seed = 0


def get_or_compute_embeddings(sequence, allow_cache=True):
    """
    Compute chunked ESM embeddings for a sequence, with optional caching.

    Returns:
        chunk_dicts: list of {"start_idx", "end_idx", "embedding" (CPU tensor)}
        full_len:    int, full sequence length
    """
    seq = sequence.strip()
    if not seq:
        return [], 0

    seq_key = seq.upper()

    if allow_cache and seq_key in embedding_cache:
        return embedding_cache[seq_key]

    chunk_dicts, length = embed_sequence_chunks(seq)

    if (
        allow_cache
        and length <= EMBEDDING_CACHE_MAX_LEN
        and len(embedding_cache) < EMBEDDING_CACHE_MAX_SEQS
    ):
        embedding_cache[seq_key] = (chunk_dicts, length)

    return chunk_dicts, length


def update_status(message, busy=False):
    print(message)


def parse_fasta(text):
    """Return list of (index, id, seq) from FASTA text."""
    records = list(SeqIO.parse(StringIO(text.strip()), "fasta"))
    return [(i, rec.id, str(rec.seq)) for i, rec in enumerate(records)]


def read_fasta(path):
    with open(path) as f:
        return f.read()


def sanitize_filename(text):
    return re.sub(r"[^a-zA-Z0-9_.-]", "_", text)


def scramble_sequence(seq, seed=None):
    rng = random.Random(seed)
    arr = list(seq)
    rng.shuffle(arr)
    return "".join(arr)


def compute_per_residue_z_and_p(real_pred, scrambled_pred_stack):
    """
    real_pred: (L,)
    scrambled_pred_stack: (N_scrambles, L)
    """
    mean_scr = np.mean(scrambled_pred_stack, axis=0)
    std_scr = np.std(scrambled_pred_stack, axis=0)
    std_scr[std_scr == 0] = 1e-6

    z_scores = (real_pred - mean_scr) / std_scr
    p_values = 1.0 - norm.cdf(z_scores)
    _, q_values = fdrcorrection(p_values, alpha=0.05)
    return z_scores, p_values, q_values


def save_predictions(
    filename,
    seq_id,
    partner_id,
    model_type,
    sequence,
    probabilities,
    p_values=None,
    q_values=None,
):
    """
    Write a tab-delimited file with per-residue probabilities and (optionally) p/FDR.
    sequence and probabilities must have the same length.
    """
    probs = np.asarray(probabilities, dtype=float)
    assert len(sequence) == len(probs), "sequence / probability length mismatch"

    with open(filename, "w") as f:
        f.write(f"# Query:{seq_id} Partner:{partner_id} Model:{model_type}-centric\n")
        if use_null_background and p_values is not None and q_values is not None:
            f.write("Position\tAminoAcid\tProbability\tP-value\tFDR\n")
            for i, (aa, prob, p, q) in enumerate(zip(sequence, probs, p_values, q_values), 1):
                f.write(f"{i}\t{aa}\t{prob:.3f}\t{p:.4g}\t{q:.4g}\n")
        else:
            f.write("Position\tAminoAcid\tProbability\n")
            for i, (aa, prob) in enumerate(zip(sequence, probs), 1):
                f.write(f"{i}\t{aa}\t{prob:.3f}\n")


def embed_sequence_chunks(sequence, max_tokens=None, stride=None):
    """
    Chunk a single amino-acid sequence, embed each chunk with ESM, and
    return a list of dicts:
      [{"start_idx": int, "end_idx": int, "embedding": [L_chunk, D]}, ...], L_full
    """
    if max_tokens is None:
        max_tokens = MAX_TOKENS
    if stride is None:
        stride = STRIDE

    esm_model.eval()
    device_local = next(esm_model.parameters()).device

    win = max_tokens - 2
    length = len(sequence)

    if length == 0:
        return [], 0

    chunk_data = []
    for start in range(0, length, stride):
        end = min(start + win, length)
        if end <= start:
            break

        enc = tokenizer(sequence[start:end], return_tensors="pt", add_special_tokens=True)
        input_ids = enc["input_ids"].to(device_local)

        with torch.no_grad():
            out = esm_model(input_ids)
            emb = out.last_hidden_state[:, 1:-1, :].squeeze(0).contiguous()

        chunk_data.append(
            {
                "start_idx": start,
                "end_idx": end,
                "embedding": emb.cpu(),
            }
        )

        if end == length:
            break

    return chunk_data, length


_ModuleBase = nn.Module if nn is not None else object


class CrossAttentionLayer(_ModuleBase):
    def __init__(self, embed_dim, num_heads):
        super().__init__()
        self.cross_attn = nn.MultiheadAttention(embed_dim, num_heads, batch_first=True)
        self.norm = nn.LayerNorm(embed_dim)

    def forward(self, query, context, context_mask=None):
        out, _ = self.cross_attn(query, context, context, key_padding_mask=context_mask)
        return self.norm(query + out)


class ChunkwiseInteractionModel(_ModuleBase):
    def __init__(self, embed_dim=None, num_heads=None, initial_bias=None):
        super().__init__()
        if embed_dim is None:
            embed_dim = EMBED_DIM
        if num_heads is None:
            num_heads = NUM_HEADS
        self.cross = CrossAttentionLayer(embed_dim, num_heads)
        self.mlp = nn.Sequential(
            nn.Linear(embed_dim, embed_dim),
            nn.ReLU(),
            nn.Dropout(0.10),
            nn.Linear(embed_dim, 1),
        )
        if initial_bias is not None:
            with torch.no_grad():
                self.mlp[-1].bias.fill_(initial_bias)

    def forward(self, q_chunks, c_chunks):
        """
        q_chunks, c_chunks: lists of [L, D] tensors (on device)
        Returns a list of [L_chunk] logit vectors, one per query chunk.
        """
        pooled_per_q = []
        for q in q_chunks:
            if q.ndim == 2:
                q = q.unsqueeze(0)
            ctx_logits = []
            for c in c_chunks:
                if c.ndim == 2:
                    c = c.unsqueeze(0)
                ctx_mask = c.abs().sum(dim=-1) == 0
                x = self.cross(q, c, context_mask=ctx_mask)
                logits = self.mlp(x).squeeze(0).squeeze(-1)
                ctx_logits.append(logits)
            pooled = torch.max(torch.stack(ctx_logits, dim=0), dim=0).values
            pooled_per_q.append(pooled)
        return pooled_per_q


def merge_chunk_logits_to_full(pooled_per_q, q_starts, full_len, device_arg=None):
    """
    pooled_per_q: list of [L_chunk] logit tensors
    q_starts: list of starting indices for each chunk
    """
    if device_arg is None:
        device_arg = device
    out = torch.full((full_len,), -1e9, dtype=torch.float32, device=device_arg)
    for vec, st in zip(pooled_per_q, q_starts):
        vec = vec.to(device_arg)
        length = vec.shape[0]
        placed = torch.full((full_len,), -1e9, dtype=torch.float32, device=device_arg)
        placed[st:st + length] = vec
        out = torch.maximum(out, placed)
    out[out <= -1e8] = -30.0
    return out


def _safe_torch_load(path, map_location=None):
    """
    Torch 2.6+ defaults to weights_only=True, which can break older checkpoints.
    This helper first tries the default behavior, and if that fails, retries with
    weights_only=False (trusted file).
    """
    try:
        return torch.load(path, map_location=map_location)
    except TypeError:
        raise
    except Exception:
        return torch.load(path, map_location=map_location, weights_only=False)


def _find_existing_model_path(primary, fallback=None):
    """
    Helper: return the first existing path among [primary, fallback].
    Raises FileNotFoundError with a clear message if nothing is found.
    """
    candidates = [primary]
    if fallback and fallback not in candidates:
        candidates.append(fallback)
    for path in candidates:
        if path and os.path.exists(path):
            return path
    raise FileNotFoundError(
        f"Some of the expected model files were not found: {candidates}. "
        "Use --download-models or provide --peptide-model-path/--receptor-model-path."
    )


def load_ppi_models():
    """
    Load peptide- and receptor-centric models using the configured model paths.
    """
    global peptide_site_model, receptor_site_model

    pep_path = _find_existing_model_path(PEPTIDE_MODEL_PATH, PEPTIDE_MODEL_FALLBACK)
    rec_path = _find_existing_model_path(RECEPTOR_MODEL_PATH, RECEPTOR_MODEL_FALLBACK)

    if peptide_site_model is None:
        print(f"Loading peptide-centric model from {pep_path}")
        m = ChunkwiseInteractionModel(embed_dim=EMBED_DIM, num_heads=NUM_HEADS).to(device)
        state = _safe_torch_load(pep_path, map_location=device)
        m.load_state_dict(state, strict=True)
        m.eval()
        peptide_site_model = m

    if receptor_site_model is None:
        print(f"Loading receptor-centric model from {rec_path}")
        m = ChunkwiseInteractionModel(embed_dim=EMBED_DIM, num_heads=NUM_HEADS).to(device)
        state = _safe_torch_load(rec_path, map_location=device)
        m.load_state_dict(state, strict=True)
        m.eval()
        receptor_site_model = m


def predict_query_given_context(query_seq, context_seq, model, cache_query=True, cache_context=True):
    """
    Predict per-residue probabilities ONLY for `query_seq`,
    using `context_seq` only as cross-attention context.

    Returns:
        probs: np.ndarray shape (len(query_seq),)
    """
    q_chunk_dicts, q_len = get_or_compute_embeddings(query_seq, allow_cache=cache_query)

    if cache_context:
        c_chunk_dicts, _ = get_or_compute_embeddings(context_seq, allow_cache=True)
    else:
        c_chunk_dicts, _ = embed_sequence_chunks(context_seq)

    q_chunks = [d["embedding"].to(device) for d in q_chunk_dicts]
    c_chunks = [d["embedding"].to(device) for d in c_chunk_dicts]
    q_starts = [int(d["start_idx"]) for d in q_chunk_dicts]

    if len(q_chunks) == 0 or len(c_chunks) == 0:
        return np.zeros(len(query_seq), dtype=float)

    with torch.no_grad():
        pooled_per_q = model(q_chunks, c_chunks)
        logits_full = merge_chunk_logits_to_full(pooled_per_q, q_starts, q_len)
        probs = torch.sigmoid(logits_full).cpu().numpy()

    if len(probs) > len(query_seq):
        probs = probs[:len(query_seq)]
    elif len(probs) < len(query_seq):
        probs = np.pad(probs, (0, len(query_seq) - len(probs)), "edge")

    return probs


def compute_scrambled_predictions_for_query(query_seq, partner_seq, model, num_scrambles=100, base_seed=0):
    """
    Null distribution by scrambling PARTNER while keeping QUERY fixed.

    - Query embeddings: cached/reused
    - Scrambled partner embeddings: NOT cached

    Returns: array shape (num_scrambles, len(query_seq))
    """
    all_scrambled = []
    for k in range(num_scrambles):
        seed = base_seed + k
        scrambled_partner = scramble_sequence(partner_seq, seed=seed)

        probs = predict_query_given_context(
            query_seq=query_seq,
            context_seq=scrambled_partner,
            model=model,
            cache_query=True,
            cache_context=False,
        )
        all_scrambled.append(probs)

    return np.stack(all_scrambled, axis=0)


def predict_query_two_views(query_seq, partner_seq, cache_query=True, cache_partner=True):
    """
    Return two probability vectors over QUERY positions only:
      - query as receptor (receptor-centric model)
      - query as ligand/peptide (peptide-centric model)

    partner_seq is context only in both cases.
    """
    p_as_receptor = predict_query_given_context(
        query_seq=query_seq,
        context_seq=partner_seq,
        model=receptor_site_model,
        cache_query=cache_query,
        cache_context=cache_partner,
    )

    p_as_peptide = predict_query_given_context(
        query_seq=query_seq,
        context_seq=partner_seq,
        model=peptide_site_model,
        cache_query=cache_query,
        cache_context=cache_partner,
    )

    return p_as_receptor, p_as_peptide


def run_ppi_for_pair(query_id, query_seq, partner_id, partner_seq, output_dir, num_scrambles=5):
    """
    For one (query, partner) pair, produce TWO probability tracks over the QUERY positions only:
      1) query-as-receptor  (receptor model)
      2) query-as-peptide   (peptide model)

    partner is context only.
    """
    print(f"  {query_id} <-> {partner_id}")

    global pred_time_total_s, pred_time_pairs, pred_num_pairs
    t0 = time.perf_counter()

    p_rec, p_pep = predict_query_two_views(query_seq, partner_seq)

    rec_pvals = rec_qvals = None
    pep_pvals = pep_qvals = None

    if use_null_background:
        print("    Computing null background (query-as-receptor)...")
        scr_rec = compute_scrambled_predictions_for_query(
            query_seq=query_seq,
            partner_seq=partner_seq,
            model=receptor_site_model,
            num_scrambles=num_scrambles,
            base_seed=base_scramble_seed,
        )
        _, rec_pvals, rec_qvals = compute_per_residue_z_and_p(p_rec, scr_rec)

        print("    Computing null background (query-as-peptide)...")
        scr_pep = compute_scrambled_predictions_for_query(
            query_seq=query_seq,
            partner_seq=partner_seq,
            model=peptide_site_model,
            num_scrambles=num_scrambles,
            base_seed=base_scramble_seed,
        )
        _, pep_pvals, pep_qvals = compute_per_residue_z_and_p(p_pep, scr_pep)

    t1 = time.perf_counter()
    dt = t1 - t0

    pred_time_total_s += dt
    pred_num_pairs += 1
    pred_time_pairs[(query_id, partner_id)] = dt

    safe_q = sanitize_filename(query_id)
    safe_p = sanitize_filename(partner_id)

    rec_file = os.path.join(output_dir, f"{safe_q}__{safe_p}__query_as_receptor.txt")
    pep_file = os.path.join(output_dir, f"{safe_q}__{safe_p}__query_as_peptide.txt")

    save_predictions(rec_file, query_id, partner_id, "query_as_receptor", query_seq, p_rec, rec_pvals, rec_qvals)
    save_predictions(pep_file, query_id, partner_id, "query_as_peptide", query_seq, p_pep, pep_pvals, pep_qvals)

    fig_name = os.path.join(output_dir, f"{safe_q}__{safe_p}__query_combined_plot.png")
    length = len(query_seq)
    x = np.arange(1, length + 1)

    fig, ax = plt.subplots(figsize=(10, 4), dpi=300)
    ax.plot(x, p_rec, label="Receptor-centric model", linewidth=2, color="orange")
    ax.plot(x, p_pep, label="Peptide-centric model", linestyle="--", linewidth=2, color="blue")

    if use_null_background and (rec_qvals is not None) and (pep_qvals is not None):
        rec_q = np.asarray(rec_qvals, dtype=float)
        pep_q = np.asarray(pep_qvals, dtype=float)

        if rec_q.shape[0] != length or pep_q.shape[0] != length:
            raise ValueError(
                f"Length mismatch: L={length}, rec_q={rec_q.shape[0]}, pep_q={pep_q.shape[0]}"
            )

        rec_sig = rec_q < 0.05
        pep_sig = pep_q < 0.05

        if np.any(rec_sig):
            ax.scatter(
                x[rec_sig],
                np.asarray(p_rec)[rec_sig],
                s=14,
                color="red",
                edgecolors="none",
                zorder=5,
                label="Significant (FDR<0.05)",
            )
        if np.any(pep_sig):
            ax.scatter(
                x[pep_sig],
                np.asarray(p_pep)[pep_sig],
                s=14,
                color="red",
                edgecolors="none",
                zorder=5,
                label=None if np.any(rec_sig) else "Significant (FDR<0.05)",
            )

    ax.set_xlabel(f"Residue position (L={length})")
    ax.set_xlim(1, length)
    ax.set_ylabel("PPI site probability")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Query: {query_id}\nPartner: {partner_id}")
    ax.grid(True, linestyle=":")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(fig_name)
    if show_plots:
        plt.show()
    plt.close(fig)


def run_ppi_predictions(query_fasta_text, partner_fasta_text, output_dir, num_scrambles=5):
    """
    Loop over all query/partner pairs from FASTA inputs and write outputs into output_dir.
    """
    seqs1 = parse_fasta(query_fasta_text)
    seqs2 = parse_fasta(partner_fasta_text)

    print(f"Parsed {len(seqs1)} query sequence(s) and {len(seqs2)} partner sequence(s).")

    if not seqs1 or not seqs2:
        print("No sequences found in one of the inputs.")
        update_status("No sequences found in one of the inputs.", busy=False)
        return

    total_pairs = len(seqs1) * len(seqs2)
    pair_idx = 0

    for (i1, query_id, query_seq), (i2, partner_id, partner_seq) in product(seqs1, seqs2):
        pair_idx += 1
        status_msg = f"Running predictions for {query_id} vs {partner_id} (pair {pair_idx}/{total_pairs})..."
        update_status(status_msg, busy=True)

        print(f"[Query {i1+1}/{len(seqs1)}] [Partner {i2+1}/{len(seqs2)}]")
        run_ppi_for_pair(
            query_id,
            query_seq,
            partner_id,
            partner_seq,
            output_dir=output_dir,
            num_scrambles=num_scrambles,
        )

        done_msg = f"Finished {query_id} vs {partner_id} (pair {pair_idx}/{total_pairs})."
        update_status(done_msg, busy=False)


def create_zip(folder, output_path=None):
    if output_path is None:
        output_path = f"{folder.rstrip(os.sep)}.zip"
    with zipfile.ZipFile(output_path, "w") as zipf:
        for file in glob(f"{folder}/*"):
            arcname = os.path.basename(file)
            zipf.write(file, arcname=arcname)
    return output_path


def collect_sites(input_dir, cutoff=0.5, fdr_limit="None"):
    current_txt_files = sorted(glob(os.path.join(input_dir, "*.txt")))
    fdr_active = fdr_limit != "None"
    output_lines = []

    if not current_txt_files:
        output_lines.append("No prediction files loaded.")
        return output_lines

    for file_path in current_txt_files:
        results = []
        try:
            with open(file_path) as f:
                header = f.readline()
                f.seek(0)
                has_fdr = "q" in header.lower() or "fdr" in header.lower()

                for line in f:
                    if line.startswith("#") or line.startswith("Position"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        continue
                    pos, aa, prob = parts[:3]
                    try:
                        prob = float(prob)
                        qval = float(parts[4]) if has_fdr and len(parts) > 4 else None
                        if prob >= cutoff and (not fdr_active or qval is None or qval < float(fdr_limit)):
                            results.append((pos, aa, prob, qval if has_fdr else None))
                    except Exception:
                        continue

            if results:
                output_lines.append(f"{os.path.basename(file_path)} - {len(results)} site(s) found:")
                header = "  Res\tProb" + ("\tFDR" if results[0][3] is not None else "")
                output_lines.append(header)
                for result in results:
                    row = f"  {result[1]}{result[0]}\t{result[2]:.3f}"
                    if result[3] is not None:
                        row += f"\t{result[3]:.3g}"
                    output_lines.append(row)
                output_lines.append("")
        except Exception as exc:
            output_lines.append(f"Error reading {file_path}: {exc}")

    return output_lines


def maybe_download_models():
    os.makedirs("models", exist_ok=True)
    if not os.path.exists(PEPTIDE_MODEL_PATH):
        print(f"Downloading peptide-centric model to {PEPTIDE_MODEL_PATH}")
        urllib.request.urlretrieve(PEPTIDE_MODEL_URL, PEPTIDE_MODEL_PATH)
    if not os.path.exists(RECEPTOR_MODEL_PATH):
        print(f"Downloading receptor-centric model to {RECEPTOR_MODEL_PATH}")
        urllib.request.urlretrieve(RECEPTOR_MODEL_URL, RECEPTOR_MODEL_PATH)


def setup_runtime(args):
    global ESM_MODEL_NAME, PEPTIDE_MODEL_URL, RECEPTOR_MODEL_URL
    global PEPTIDE_MODEL_PATH, RECEPTOR_MODEL_PATH
    global NUM_HEADS, MAX_TOKENS, STRIDE, PRED_CUTOFF
    global EMBEDDING_CACHE_MAX_SEQS, EMBEDDING_CACHE_MAX_LEN
    global device, tokenizer, esm_model, EMBED_DIM
    global use_null_background, show_plots, base_scramble_seed

    ESM_MODEL_NAME = args.esm_model_name
    PEPTIDE_MODEL_URL = args.peptide_model_url
    RECEPTOR_MODEL_URL = args.receptor_model_url
    PEPTIDE_MODEL_PATH = args.peptide_model_path
    RECEPTOR_MODEL_PATH = args.receptor_model_path
    NUM_HEADS = args.num_heads
    MAX_TOKENS = args.max_tokens
    STRIDE = args.stride
    PRED_CUTOFF = args.prediction_cutoff
    EMBEDDING_CACHE_MAX_SEQS = args.embedding_cache_max_seqs
    EMBEDDING_CACHE_MAX_LEN = args.embedding_cache_max_len
    use_null_background = args.use_null_background
    show_plots = args.show_plots
    base_scramble_seed = args.seed

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    print("Using device:", device, "\n")

    if args.download_models:
        maybe_download_models()

    print(f"Loading ESM-2 checkpoint: {ESM_MODEL_NAME}...")
    tokenizer = EsmTokenizer.from_pretrained(ESM_MODEL_NAME)
    esm_model = EsmModel.from_pretrained(ESM_MODEL_NAME).to(device)

    for param in esm_model.parameters():
        param.requires_grad = False
    esm_model.eval()
    EMBED_DIM = esm_model.config.hidden_size

    print("\nESM-2 checkpoint loaded")
    print(
        "Note: You may safely ignore the warning about uninitialized pooler weights above. "
        "This component is not used in the prediction workflow and does not affect results.\n"
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run SPPIDER-seq partner-aware PPI site predictions from FASTA files."
    )
    parser.add_argument("--query", required=True, help="FASTA file containing query sequence(s).")
    parser.add_argument("--partners", required=True, help="FASTA file containing partner sequence(s).")
    parser.add_argument("--output-dir", required=True, help="Output directory for prediction files.")

    parser.add_argument("--use-null-background", action="store_true", help="Estimate significance with scrambled partners.")
    parser.add_argument("--num-scrambles", type=int, default=5, help="Scrambles for null background. Default: 5")
    parser.add_argument("--show-plots", action="store_true", help="Show plots as they are generated.")
    parser.add_argument("--min-probability", type=float, default=0.5, help="Minimum probability for site summary. Default: 0.5")
    parser.add_argument(
        "--fdr-cutoff",
        default="None",
        choices=["None", "0.05", "0.01", "0.001"],
        help="FDR cutoff for site summary. Default: None",
    )

    parser.add_argument("--esm-model-name", default=ESM_MODEL_NAME, help=f"ESM model name. Default: {ESM_MODEL_NAME}")
    parser.add_argument("--peptide-model-path", default=PEPTIDE_MODEL_PATH, help=f"Peptide model path. Default: {PEPTIDE_MODEL_PATH}")
    parser.add_argument("--receptor-model-path", default=RECEPTOR_MODEL_PATH, help=f"Receptor model path. Default: {RECEPTOR_MODEL_PATH}")
    parser.add_argument("--peptide-model-url", default=PEPTIDE_MODEL_URL, help="Peptide model download URL.")
    parser.add_argument("--receptor-model-url", default=RECEPTOR_MODEL_URL, help="Receptor model download URL.")
    parser.add_argument("--download-models", action="store_true", help="Download missing notebook-style model files before prediction.")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto", help="Inference device. Default: auto")
    parser.add_argument("--num-heads", type=int, default=16, help="Cross-attention heads. Default: 16")
    parser.add_argument("--max-tokens", type=int, default=1024, help="ESM token chunk size including BOS/EOS. Default: 1024")
    parser.add_argument("--stride", type=int, default=512, help="Residue stride for chunked embeddings. Default: 512")
    parser.add_argument("--prediction-cutoff", type=float, default=0.50, help="Plotting/site cutoff. Default: 0.50")
    parser.add_argument("--embedding-cache-max-seqs", type=int, default=256, help="Maximum cached sequence count. Default: 256")
    parser.add_argument("--embedding-cache-max-len", type=int, default=3000, help="Maximum sequence length eligible for caching. Default: 3000")
    parser.add_argument("--seed", type=int, default=0, help="Base seed for scrambled null-background generation. Default: 0")

    parser.add_argument("--make-zip", action="store_true", help="Create a ZIP archive for the output directory.")
    parser.add_argument("--sites-output", help="Optional path for a text summary of residues passing filters.")
    return parser


def main():
    args = build_parser().parse_args()

    if IMPORT_ERROR is not None:
        raise SystemExit(
            f"Missing required dependency: {IMPORT_ERROR.name}. "
            "Install the notebook requirements before running predictions."
        )

    os.makedirs(args.output_dir, exist_ok=True)

    global pred_time_total_s, pred_time_pairs, pred_num_pairs
    pred_time_total_s = 0.0
    pred_time_pairs = {}
    pred_num_pairs = 0
    embedding_cache.clear()

    setup_runtime(args)
    load_ppi_models()

    query_fasta_text = read_fasta(args.query)
    partner_fasta_text = read_fasta(args.partners)

    print(f"Output folder: {args.output_dir}")
    run_ppi_predictions(
        query_fasta_text=query_fasta_text,
        partner_fasta_text=partner_fasta_text,
        output_dir=args.output_dir,
        num_scrambles=args.num_scrambles,
    )

    if pred_num_pairs > 0:
        avg = pred_time_total_s / pred_num_pairs
        msg = (
            f"All predictions completed. "
            f"Prediction compute time: {pred_time_total_s:.2f}s total "
            f"({avg:.2f}s per pair; {pred_num_pairs} pair(s))."
        )
    else:
        msg = "All predictions completed."
    update_status(msg, busy=False)

    if args.sites_output:
        site_lines = collect_sites(args.output_dir, args.min_probability, args.fdr_cutoff)
        with open(args.sites_output, "w") as f:
            f.write("\n".join(site_lines))
        print(f"Wrote {args.sites_output}")

    if args.make_zip:
        zip_path = create_zip(args.output_dir)
        print(f"Wrote {zip_path}")


if __name__ == "__main__":
    main()
