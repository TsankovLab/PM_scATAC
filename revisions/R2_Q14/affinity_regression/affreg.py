"""
Affinity regression (Pelossof et al. 2015; Osmanbeyoglu et al. 2017), ported line-by-line
from the Osmanbeyoglu lab MATLAB code (github.com/osmanbeyoglulab/Affreg):
  ar_train.m, ar_model2w.m, ar_predict.m, ar_reconstruction.m, SLEP/LeastR.m (+initFactor, sll_opts)

Model:  Y ~ D W P'
  Y  genes x samples        (mean-centred expression)
  D  genes x TFs            (binary TF-target prior)
  P  samples x proteins     (mean-centred RPPA)
  W  TFs x proteins         (learned interaction matrix)
TF activity per sample      = W P'          (TFs x samples)
(phospho)protein activity   = Y' D W        (samples x proteins)

MATLAB is column-major: every vec()/reshape()/linear index below uses order='F' to match.
Ported because the cluster MATLAB licence excludes this user; the port is validated on
simulated data with a known W (run this file directly).
"""
import numpy as np
from scipy.linalg import orth


def least_r(A, y, z, rsL2=0.0, tol=1e-3, max_iter=10000, plus001=True):
    """SLEP LeastR with sll_opts defaults (init=0, tFlag=0, nFlag=0, rFlag=0, mFlag=0, lFlag=0):
       min 1/2||Ax-y||^2 + 1/2 rsL2 ||x||^2 + z ||x||_1
       Accelerated gradient with Armijo-Goldstein line search.  plus001 mirrors the
       '+0.001' present in the lab's copy of LeastR.m."""
    m, n = A.shape
    ATy = A.T @ y
    lam = z
    x = ATy.copy()                     # init 0: x0 = ATy, rescaled by initFactor
    Ax = A @ x
    x_norm, x_2norm = np.abs(x).sum(), x @ x
    if x_norm >= 1e-6:
        ratio = (Ax @ y - lam * x_norm) / (Ax @ Ax + rsL2 * x_2norm)
        x = ratio * x
        Ax = ratio * Ax
    L = 1.0 + rsL2
    xp, Axp, xxp = x.copy(), Ax.copy(), np.zeros(n)
    alphap, alpha = 0.0, 1.0
    fun = []
    for it in range(max_iter):
        beta = (alphap - 1.0) / alpha
        s = x + beta * xxp
        As = Ax + beta * (Ax - Axp)
        g = A.T @ As - ATy + rsL2 * s
        xp, Axp = x, Ax
        bflag = False
        while True:
            v = s - g / L
            x = np.sign(v) * np.maximum(np.abs(v) - lam / L, 0.0)
            v = x - s
            Ax = A @ x
            Av = Ax - As
            r_sum, l_sum = v @ v, Av @ Av
            if r_sum <= 1e-20:
                bflag = True
                break
            if l_sum <= r_sum * (L - rsL2) + (0.001 if plus001 else 0.0):
                break
            L = max(2.0 * L, l_sum / r_sum + rsL2)
        alphap, alpha = alpha, (1.0 + np.sqrt(4.0 * alpha * alpha + 1.0)) / 2.0
        xxp = x - xp
        Axy = Ax - y
        fun.append(Axy @ Axy / 2.0 + rsL2 / 2.0 * (x @ x) + np.abs(x).sum() * lam)
        if bflag:
            break
        if it >= 1 and abs(fun[-1] - fun[-2]) <= tol:
            break
    return x, np.array(fun)


def _first_ge(cum, frac):
    ix = np.nonzero(cum >= frac - 1e-12)[0]
    return int(ix[0]) + 1


def ar_train(D, P, Y, lam, rsL2=0.0, spectrumA=1.0, spectrumB=1.0, old_version=False,
             tol=1e-3, max_iter=10000):
    """ar_train.m.  D genes x TFs, P samples x proteins, Y genes x samples."""
    A = Y.T @ D                         # samples x TFs
    B = P.T                             # proteins x samples
    Z = Y.T @ Y                         # samples x samples
    UA, sA, VAt = np.linalg.svd(A, full_matrices=True)
    UB, sB, VBt = np.linalg.svd(B, full_matrices=True)
    da = _first_ge(np.cumsum(sA) / sA.sum(), spectrumA)
    db = _first_ge(np.cumsum(sB) / sB.sum(), spectrumB)
    Ua, Va = UA[:, :da], VAt.T[:, :da]
    Ub, Vb = UB[:, :db], VBt.T[:, :db]
    L = np.kron(Vb, Ua)                 # (M*M) x (db*da)
    Yv = Z.flatten(order="F")
    if not old_version:                 # drop the diagonal (self-similarity) equations
        M = Z.shape[0]
        diag_ix = np.arange(M) * M + np.arange(M)
        keep = np.ones(M * M, dtype=bool); keep[diag_ix] = False
        L, Yv = L[keep], Yv[keep]
    beta, fun = least_r(L, Yv, lam, rsL2, tol=tol, max_iter=max_iter)
    return dict(Ua=Ua, Ub=Ub, Sa=sA[:da], Sb=sB[:db], Va=Va, Vb=Vb, beta=beta,
                lam=lam, rsL2=rsL2, da=da, db=db, n_iter=len(fun), fun=fun)


def _pinv_diag(s):
    tol = max(s.shape) * np.finfo(float).eps * (s.max() if s.size else 0)
    return np.diag(np.where(s > tol, 1.0 / s, 0.0))


def ar_model2w(model):
    """ar_model2w.m: W = Va pinv(Sa) reshape(beta, da, db) pinv(Sb) Ub'"""
    X = model["beta"].reshape((model["da"], model["db"]), order="F")
    return model["Va"] @ _pinv_diag(model["Sa"]) @ X @ _pinv_diag(model["Sb"]) @ model["Ub"].T


def ar_reconstruction(Y_train, pred_test):
    A = Y_train.T @ pred_test
    O = orth(Y_train)
    c = np.linalg.lstsq(O, Y_train, rcond=None)[0]
    ct = np.linalg.lstsq(c.T, A, rcond=None)[0]
    return O @ ct


def ar_predict(D, P_test, Y_train, model):
    W = ar_model2w(model)
    pred = D @ (W @ P_test.T)
    return dict(rec=ar_reconstruction(Y_train, pred), pred=pred, W=W)


def tf_activity(W, P):      return W @ P.T          # TFs x samples
def protein_activity(Y, D, W): return Y.T @ D @ W   # samples x proteins


if __name__ == "__main__":
    # mirror run_helloworld.m at reduced size: noiseless Y = D W P', known W
    rng = np.random.default_rng(1)
    n, m, p, q = 300, 70, 40, 30
    D, P, W = rng.standard_normal((n, p)), rng.standard_normal((m, q)), rng.standard_normal((p, q))
    idx = rng.permutation(m); tr, te = idx[:49], idx[49:]
    Ytr, Yte = D @ W @ P[tr].T, D @ W @ P[te].T
    for lam, sA, sB in [(0.1, 0.95, 0.9), (0.01, 1.0, 1.0), (0.001, 1.0, 1.0)]:
        mod = ar_train(D, P[tr], Ytr, lam, 0.0, sA, sB)
        pr = ar_predict(D, P[te], Ytr, mod)
        cw = np.corrcoef(W.ravel(), pr["W"].ravel())[0, 1]
        cy = np.mean([np.corrcoef(Yte[:, j], pr["rec"][:, j])[0, 1] for j in range(len(te))])
        ca = np.corrcoef(tf_activity(W, P[te]).ravel(), tf_activity(pr["W"], P[te]).ravel())[0, 1]
        print(f"lambda={lam:<6} specA={sA} specB={sB}: da={mod['da']} db={mod['db']} iters={mod['n_iter']} | "
              f"corr(W)={cw:.4f} | mean corr(Y_test, rec)={cy:.4f} | corr(TF activity)={ca:.4f}")
