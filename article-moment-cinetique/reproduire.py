"""
Reproduit l'ensemble des nombres et des figures de l'article
« La Lune et l'évacuation du moment cinétique » (versions FR et EN).

    python3 reproduire.py            -> resultats.json + 6 figures PNG
    python3 reproduire.py --figures  -> redessine les figures sans recalculer

Dépend de stripping.py (moteur) et analyse.py (classement), numpy, scipy,
matplotlib. Durée : quelques minutes sur un portable.

M. Debailleul, 2026 -- CC BY-NC 4.0
"""
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from stripping import Config, ME, ML, RE, LEM, GM, AROCHE
from analyse import run_case

# --------------------------------------------------------------- simulations
REF = {
    "coplanaire i=0":  dict(iB_deg=0.0),
    "prograde i=20":   dict(iB_deg=20.0),
    "prograde i=25":   dict(iB_deg=25.0),
    "prograde i=30":   dict(iB_deg=30.0),
    "retrograde i=25": dict(iB_deg=25.0, retrograde=True),
}
CATS = ["fenetre_lunaire", "disque_interne", "retombee", "echappee", "emportee_B"]


def released_j(r):
    """moment cinetique specifique (vecteur) de chaque element au lacher,
    dans le repere inertiel : j = Omega (varpi^2 z_hat - z r_perp)."""
    cfg = r["cfg"]
    ang = cfg.Omega * np.nan_to_num(r["t_rel"])
    c, s = np.cos(ang), np.sin(ang)
    p0 = r["pos0"]
    x = c * p0[:, 0] - s * p0[:, 1]
    y = s * p0[:, 0] + c * p0[:, 1]
    z = p0[:, 2]
    return cfg.Omega * np.stack([-z * x, -z * y, x * x + y * y], axis=1)


def band_stat(r, q):
    b = np.isin(r["cat"], ["fenetre_lunaire", "disque_interne"]) & \
        (np.abs(r["lat0"]) < np.radians(5))
    return float(np.sqrt(np.nanpercentile(r["rcirc"][b], q))) if b.any() else None


def simulate():
  print("Cas de reference (grille 41 x 240)...", flush=True)
  res = {n: run_case(n, **kw) for n, kw in REF.items()}

  out = {"reference": {}, "scan": {}}
  for n, r in res.items():
      w, cat = r["wmass"], r["cat"]
      lifted = cat != "non_leve"
      L = w[lifted].sum()
      d = {"souleve_bande": float(L)}
      d.update({c: float(w[cat == c].sum() / L) for c in CATS})
      d["L_d9"] = band_stat(r, 90)
      d["L_med"] = band_stat(r, 50)

      cfg = r["cfg"]
      unit = cfg.Omega * cfg.Req**2
      j = released_j(r)
      cap = cat == "emportee_B"
      if cap.any():
          J = (j[cap] * w[cap, None]).sum(0)
          d["jz_capture"] = float(J[2] / w[cap].sum() / unit)
          d["angle_vecteur_emporte_deg"] = float(np.degrees(np.arctan2(np.hypot(J[0], J[1]), J[2])))
      d["jz_souleve"] = float((j[lifted, 2] * w[lifted]).sum() / L / unit)
      t = r["t_rel"][r["released"]] / 3600.0
      d["lachers_dans_pm2h"] = float(np.mean(np.abs(t - np.median(t)) < 2.0))
      d["etendue_lachers_h"] = float(t.max() - t.min())
      out["reference"][n] = d
      print(f"  {n:16s}", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in d.items()}, flush=True)

  print("Balayage de parametres (grille 21 x 120)...", flush=True)
  for a in [1.70, 1.80]:
      for mu in [0.50, 0.75, 0.95]:
          for b in [2.8, 3.0, 3.5]:
              if b < a + mu**0.27 + 0.02:          # contact
                  continue
              for i in [20.0, 25.0, 30.0]:
                  if i != 25.0 and (mu != 0.95 or a != 1.70):
                      continue
                  r = run_case("x", nl=21, nn=120, a=a, muB=mu, b=b, iB_deg=i)
                  w, cat = r["wmass"], r["cat"]
                  L = w[cat != "non_leve"].sum()
                  cap = cat == "emportee_B"
                  cfg = r["cfg"]
                  j = released_j(r)
                  jz = float((j[cap, 2] * w[cap]).sum() / w[cap].sum() / (cfg.Omega * cfg.Req**2)) if cap.any() else None
                  lif = cat != "non_leve"
                  jzl = float((j[lif, 2] * w[lif]).sum() / L / (cfg.Omega * cfg.Req**2))
                  key = f"a={a:.2f} mu={mu:.2f} q={b:.1f} i={i:.0f}"
                  out["scan"][key] = dict(souleve=float(L), phi=float(w[cap].sum() / L), jz=jz,
                                          jz_souleve=jzl,
                                          psi=float(w[cap].sum() / L) * (jz / jzl if jz else 0.0))
                  print(f"  {key}  souleve={L:.3f}  phi={w[cap].sum()/L:.3f}  jz={jz}", flush=True)

  # ------------------------------------------------------------ bilan global
  cfg = Config()
  S = 0.19 * ME * cfg.Omega * cfg.Req**2 / LEM
  phis = [v["phi"] for k, v in out["scan"].items() if v["souleve"] > 0.05]
  jzs = [v["jz"] for v in out["scan"].values() if v["jz"] is not None]
  phi_ref = out["reference"]["prograde i=25"]["emportee_B"]
  jz_ref = out["reference"]["prograde i=25"]["jz_capture"]
  jzl_ref = out["reference"]["prograde i=25"]["jz_souleve"]
  ratio = jz_ref / jzl_ref
  psis = [v["psi"] for v in out["scan"].values() if v["souleve"] > 0.05]
  unit = cfg.Omega * cfg.Req**2
  out["bilan"] = dict(
      spin_LEM=S,
      exces_LEM=S - 1.0,
      Omega_Req2=unit,
      jz_ref_SI=jz_ref * unit,
      masse_pour_fermeture_ML=(S - 1.0) * LEM / (jz_ref * unit) / ML,
      jz_souleve_ref=jzl_ref,
      rapport_jz_capture_sur_souleve=ratio,
      borne_naive_phiS_LEM=phi_ref * S,
      borne_retrait_LEM=phi_ref * ratio * S,
      borne_fraction_du_besoin=phi_ref * ratio * S / (S - 1.0),
      borne_retrait_grille_max_LEM=max(psis) * S,
      masse_soulevable_max_ML=0.19 * ME / jzl_ref / ML,
      masse_capturee_max_ML=phi_ref * 0.19 * ME / jzl_ref / ML,
      phi_min=min(phis), phi_max=max(phis),
      jz_min=min(jzs), jz_max=max(jzs),
  )
  for m in (30.0, 42.0):
      out["bilan"][f"moment_de_{int(m)}ML_LEM"] = m * ML * jzl_ref * unit / LEM
  print("Bilan :", {k: round(v, 4) for k, v in out["bilan"].items()}, flush=True)
  json.dump(out, open("resultats.json", "w"), indent=1)
  r = res["prograde i=25"]; cap = r["cat"] == "emportee_B"
  np.savez("histogramme_jz.npz",
           jz=released_j(r)[cap, 2] / (r["cfg"].Omega * r["cfg"].Req**2),
           w=r["wmass"][cap])
  return out


if "--figures" in sys.argv:
    out = json.load(open("resultats.json"))
else:
    out = simulate()
H = np.load("histogramme_jz.npz")

# ------------------------------------------------------------------ figures
# Palette : ordre categoriel fixe, valide daltonisme (voir note en fin de fichier)
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GRAY, INK, INK2 = "#b4b2ab", "#0b0b0b", "#52514e"

plt.rcParams.update({"font.size": 9, "figure.dpi": 160, "axes.grid": True,
                     "grid.alpha": .25, "grid.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": INK2, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})

T = {
 "en": dict(dec=".", q=r"closest approach $q_B$  ($R_\oplus$)",
    phi=r"$\phi_{\oplus\to B}$  (captured / lifted)", t1="Captured fraction stays near one quarter",
    gate="viability gate 0.20", jz=r"$j_z/\Omega_\oplus R_{\rm eq}^2$ at release",
    mass="mass (arbitrary units)", t2="Specific angular momentum removed",
    inc=r"passage inclination $i_B$ (deg)", rc=r"$\sqrt{r_{\rm circ}/R_\oplus}$",
    t3="Lifted matter does not reach the Roche limit", cos=r"closed form with $\cos^2 i$",
    roche="Roche threshold 1.704", rot="rotation alone 1.265",
    d9="simulation, upper decile", med="simulation, median", win="geometric\nwindow",
    t4="Fate of the lifted matter",
    cats=["captured by B", r"beyond $a_{\rm R}$", "inner disk", "escaped", "fell back"],
    xt=[r"$i=0^\circ$", r"$20^\circ$", r"$25^\circ$", r"$30^\circ$", "retro.\n" + r"$25^\circ$"],
    frac="fraction of lifted matter", t5="Angular-momentum budget",
    bud=["excess to\nremove", "upper bound for\nthis channel", "tidal torque on\nthe figure"],
    ylab=r"angular momentum  ($L_{\rm EM}$)"),
 "fr": dict(dec=",", q=r"distance au périastre $q_B$  ($R_\oplus$)",
    phi=r"$\phi_{\oplus\to B}$  (capturée / soulevée)", t1="La fraction capturée reste proche d'un quart",
    gate="porte de viabilité 0,20", jz=r"$j_z/\Omega_\oplus R_{\rm eq}^2$ au lâcher",
    mass="masse (unités arbitraires)", t2="Moment cinétique spécifique emporté",
    inc=r"inclinaison du passage $i_B$ (degrés)", rc=r"$\sqrt{r_{\rm circ}/R_\oplus}$",
    t3="La matière soulevée n'atteint pas la limite de Roche", cos=r"forme fermée en $\cos^2 i$",
    roche="seuil de Roche 1,704", rot="rotation seule 1,265",
    d9="simulation, décile supérieur", med="simulation, médiane", win="fenêtre\ngéométrique",
    t4="Devenir de la matière soulevée",
    cats=["capturée par B", r"au-delà de $a_{\rm R}$", "disque interne", "échappée", "retombée"],
    xt=[r"$i=0^\circ$", r"$20^\circ$", r"$25^\circ$", r"$30^\circ$", "rétro.\n" + r"$25^\circ$"],
    frac="fraction de la matière soulevée", t5="Bilan de moment cinétique",
    bud=["excédent à\névacuer", "borne supérieure\nde ce canal", "couple de marée\nsur la figure"],
    ylab=r"moment cinétique  ($L_{\rm EM}$)"),
}


def num(x, lg, nd=2, math=False):
    s = f"{x:.{nd}f}"
    if lg != "fr":
        return s
    return s.replace(".", "{,}") if math else s.replace(".", ",")


from matplotlib.ticker import FuncFormatter


def fr_ticks(ax, lg, which="y", nd=None):
    if lg != "fr":
        return
    def f(v, _):
        t = (f"{v:.{nd}f}" if nd is not None else f"{v:g}")
        return t.replace(".", ",")
    if "y" in which:
        ax.yaxis.set_major_formatter(FuncFormatter(f))
    if "x" in which:
        ax.xaxis.set_major_formatter(FuncFormatter(f))


K = np.sqrt(2) * 0.95 / np.sqrt(1.95)
ROT, TID = 0.97 * np.sqrt(1.70), K * 1.70**2 / 3.0**1.5
names = ["coplanaire i=0", "prograde i=20", "prograde i=25", "prograde i=30"]
inc = np.array([0, 20, 25, 30])

for lg, t in T.items():
    # ---- figure 1 : robustesse de la capture + moment emporte
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.5))
    qs = [2.8, 3.0, 3.5]
    for a, ls, mk in [(1.70, "-", "o"), (1.80, "--", "s")]:
        for mu, c in [(0.50, BLUE), (0.75, ORANGE), (0.95, AQUA)]:
            y = [out["scan"].get(f"a={a:.2f} mu={mu:.2f} q={q:.1f} i=25", {}).get("phi", np.nan) for q in qs]
            ax[0].plot(qs, y, ls, marker=mk, color=c, ms=5, lw=1.4,
                       label=rf"$\mu_B={num(mu, lg, math=True)}$, $a={num(a, lg, math=True)}$")
    ax[0].axhspan(0.21, 0.28, color=INK2, alpha=.08, lw=0)
    ax[0].axhline(0.20, color=INK, ls=":", lw=1.1)
    ax[0].text(2.82, 0.193, t["gate"], fontsize=7, ha="left", va="top", color=INK2)
    ax[0].set_xlabel(t["q"]); ax[0].set_ylabel(t["phi"])
    ax[0].set_ylim(0, 0.34); ax[0].set_xticks(qs)
    ax[0].set_xticklabels([num(q, lg, 1) for q in qs])
    ax[0].legend(fontsize=6.4, ncol=2, loc="lower left")
    ax[0].set_title(t["t1"], fontsize=9.5)
    fr_ticks(ax[0], lg, "y", 2)

    jz, wj = H["jz"], H["w"]
    ax[1].hist(jz, bins=40, weights=wj, color=BLUE, edgecolor="white", lw=0.6)
    m = np.average(jz, weights=wj)
    ax[1].axvline(m, color=INK, lw=1.2)
    ax[1].text(m - 0.004, ax[1].get_ylim()[1] * .92,
               r"$\langle j_z\rangle=" + num(m, lg, math=True) + r"$", ha="right", fontsize=8)
    ax[1].set_xlabel(t["jz"]); ax[1].set_ylabel(t["mass"])
    ax[1].set_title(t["t2"], fontsize=9.5)
    fr_ticks(ax[1], lg, "xy")
    plt.tight_layout(); plt.savefig(f"fig_capture_{lg}.png", bbox_inches="tight"); plt.close()

    # ---- figure 2 : controle negatif + devenir de la matiere
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.5))
    ii = np.linspace(0, 35, 200)
    ax[0].fill_between([18, 27], 1.2, 2.0, color=INK2, alpha=.08, lw=0)
    ax[0].text(22.5, 1.94, t["win"], ha="center", va="top", fontsize=7, color=INK2)
    ax[0].plot(ii, ROT + TID * np.cos(np.radians(ii))**2, color=BLUE, lw=1.6, label=t["cos"])
    ax[0].axhline(1.704, color=INK, ls="--", lw=1.1, label=t["roche"])
    ax[0].axhline(ROT, color=INK2, ls="-.", lw=1.0, label=t["rot"])
    ax[0].plot(inc, [out["reference"][n]["L_d9"] for n in names], "o-", color=ORANGE, ms=5, lw=1.4, label=t["d9"])
    ax[0].plot(inc, [out["reference"][n]["L_med"] for n in names], "s-", color=AQUA, ms=5, lw=1.4, label=t["med"])
    ax[0].set_xlabel(t["inc"]); ax[0].set_ylabel(t["rc"]); ax[0].set_ylim(1.2, 1.98)
    ax[0].legend(fontsize=6.4, loc="lower left"); ax[0].set_title(t["t3"], fontsize=9.5)
    fr_ticks(ax[0], lg, "y", 1)

    order = ["emportee_B", "fenetre_lunaire", "disque_interne", "echappee", "retombee"]
    cols = [BLUE, ORANGE, AQUA, YELLOW, GRAY]
    names2 = names + ["retrograde i=25"]
    bot = np.zeros(len(names2))
    for c, l, co in zip(order, t["cats"], cols):
        v = np.array([out["reference"][n][c] for n in names2])
        ax[1].bar(range(len(names2)), v, bottom=bot, label=l, color=co, width=.62,
                  edgecolor="white", lw=1.0)
        bot += v
    ax[1].set_xticks(range(len(names2))); ax[1].set_xticklabels(t["xt"], fontsize=8)
    ax[1].set_ylabel(t["frac"]); ax[1].set_ylim(0, 1.0)
    h, l = ax[1].get_legend_handles_labels()
    ax[1].legend(h[::-1], l[::-1], fontsize=7, loc="center left", bbox_to_anchor=(1.01, .5))
    ax[1].set_title(t["t4"], fontsize=9.5)
    ax[1].grid(axis="x", visible=False)
    fr_ticks(ax[1], lg, "y", 1)
    plt.tight_layout(); plt.savefig(f"fig_controle_{lg}.png", bbox_inches="tight"); plt.close()

    # ---- figure 3 : bilan
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    vals = [out["bilan"]["exces_LEM"], out["bilan"]["borne_retrait_LEM"], 0.07]
    ax.bar(range(3), vals, color=[INK2, BLUE, GRAY], width=.55)
    ax.errorbar(2, 0.07, yerr=[[0.05], [0.05]], color=INK, capsize=4, lw=1)
    tops = [vals[0], vals[1], 0.12]
    for k, (v, top) in enumerate(zip(vals, tops)):
        ax.text(k, top + 0.03, num(v, lg), ha="center", fontsize=9)
    ax.set_xticks(range(3)); ax.set_xticklabels(t["bud"], fontsize=7.5)
    ax.set_ylabel(t["ylab"]); ax.set_ylim(0, 1.25); ax.set_title(t["t5"], fontsize=9.5)
    ax.grid(axis="x", visible=False)
    fr_ticks(ax, lg, "y", 1)
    plt.tight_layout(); plt.savefig(f"fig_bilan_{lg}.png", bbox_inches="tight"); plt.close()

print("Figures ecrites.")

# Note sur la palette : bleu, orange, aqua, jaune dans cet ordre (palette de
# reference validee pour la separation daltonienne des paires adjacentes),
# gris neutre pour la categorie « retombee » ; chaque couleur est doublee
# d'une legende ou d'une etiquette d'axe.
