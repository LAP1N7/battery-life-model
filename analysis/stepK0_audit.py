"""K0: KIT 데이터 점검 (트랙 B, 과업 B2) -> out/K/K0_audit.md 를 자동으로 채움

[무엇을 하나]  KIT 결과 데이터(EOC, EIS)를 훑어서 사전 등록 시험을 할 수 있는지 4가지를 점검한다.
  ① SOH 0.85 도달 셀 수 (10개 이상이면 진행)   ② check-up별 시간·누적 충전량(EFC 계산용) 열이 있는가
  ③ 첫 check-up에 상온(약 25 °C) EIS가 있는가   ④ 셀별 조건표 (파일 이름의 파라미터 ID + cfg 파일)

[확인된 사실 — 데이터 논문 (Luh & Blank, Sci Data 2024) 기준]
  - 파일 이름: cell_eocv2_P<ID>_<nr>_S<slave>_C<ch>.csv (EOC), cell_eisv2_... (EIS), cell_plsv2_... (PULSE), cell_cfg_... (설정)
  - EOC 열 중 이름이 확인된 것: cap_aged_est_Ah (남은 용량), soh_cap (공칭 용량 기준 SOH), timestamp (Unix 초)
  - EIS 열 중 이름이 확인된 것: z_re_comp_mOhm, z_im_comp_mOhm, z_amp_mOhm, z_ph_deg (+ 주파수, 공칭 SoC, 조건 열)
[확인 안 된 것 — 그래서 이 스크립트가 자동으로 찾고, 못 찾으면 "확인 필요"로 적음]
  - 구분자(쉼표/세미콜론), 나머지 열의 정확한 이름, check-up(CU)·상온(RT)을 나타내는 조건 열의 값
  -> 첫 실행은 "탐색"이다. 출력된 열 이름 목록(K0_columns.txt)을 A에게 공유하면 열 이름을 고정한다.

[SOH 정의 — 사전 등록과 동일]  각 check-up 용량 / 첫 check-up 용량  (논문의 soh_cap은 공칭 기준이라 참고로만 표시)

입력: Mothat/KIT/ 아래 압축을 푼 CSV (하위 폴더 포함 자동 탐색)
출력: out/K/K0_audit.md (템플릿을 자동 채움), out/K/K0_conditions.csv, out/K/K0_columns.txt
변경 이력: 2026-11-19 최초 작성 (Claude). 실제 KIT 파일로는 실행해 보지 못함 -> 첫 실행 출력으로 검증 필요
"""
import re, glob, os
import numpy as np, pandas as pd
from common import OUT

KIT = OUT.parent.parent / "KIT"          # Mothat/KIT
KOUT = OUT / "K"; KOUT.mkdir(exist_ok=True)
SOH_TH = 0.85
RT_RANGE = (20.0, 30.0)                   # 상온 판정: EIS 평균 온도가 이 범위 [°C]

# ---------- 공통 도우미 ----------
def find_files(pattern):
    return sorted(glob.glob(str(KIT / "**" / pattern), recursive=True))

def read_csv_any(path, nrows=None):
    """구분자를 자동으로 추정해서 읽음 (쉼표/세미콜론/탭)"""
    return pd.read_csv(path, sep=None, engine="python", nrows=nrows)

def pick(cols, *keys, avoid=()):
    """열 이름에 keys가 모두 들어 있는 첫 열 (대소문자 무시). 없으면 None"""
    for c in cols:
        lc = c.lower()
        if all(k in lc for k in keys) and not any(a in lc for a in avoid): return c
    return None

def cell_id(path):
    m = re.search(r"P(\d{3})_(\d+)", os.path.basename(path))
    return (f"P{m.group(1)}_{m.group(2)}", int(m.group(1))) if m else (os.path.basename(path), None)

lines_cols, notes = [], []

# ---------- 0. 파일 목록 ----------
eoc_files, eis_files, cfg_files = find_files("cell_eoc*.csv"), find_files("cell_eis*.csv"), find_files("cell_cfg*.csv")
archives = find_files("*.zip") + find_files("*.7z")
print(f"[0] KIT 폴더: {KIT}")
print(f"    EOC {len(eoc_files)}개, EIS {len(eis_files)}개, CFG {len(cfg_files)}개, 압축 파일 {len(archives)}개")
if not eoc_files:
    print("    EOC CSV가 없음 -> KIT 폴더에 압축을 풀었는지 확인 (zip/7z는 직접 풀어야 함). 여기서 멈춤."); raise SystemExit
if archives: notes.append(f"압축 파일 {len(archives)}개가 남아 있음 (풀었는지 확인)")

for kind, files in [("EOC", eoc_files), ("EIS", eis_files), ("CFG", cfg_files)]:
    if not files: continue
    d = read_csv_any(files[0], nrows=5)
    lines_cols.append(f"=== {kind}: {os.path.basename(files[0])} ({d.shape[1]}열) ===")
    lines_cols += [f"  {c}  [{d[c].dtype}]  예: {d[c].iloc[0] if len(d) else ''}" for c in d.columns]
(KOUT / "K0_columns.txt").write_text("\n".join(lines_cols), encoding="utf-8")
print(f"    열 이름 목록 저장: {KOUT / 'K0_columns.txt'}")

# ---------- ①② EOC: 용량·시간·누적 충전량 ----------
d0 = read_csv_any(eoc_files[0], nrows=5); cols = list(d0.columns)
c_cap = pick(cols, "cap_aged_est") or pick(cols, "cap", "ah")
c_soh = pick(cols, "soh_cap")
c_ts = pick(cols, "timestamp", avoid=("origin",))
c_cond = pick(cols, "cyc_cond") or pick(cols, "condition")
c_q = pick(cols, "total", "q", avoid=("cu", "rt", "ot", "e_")) or pick(cols, "q_total") or pick(cols, "total_q")
print(f"[①②] 자동으로 찾은 EOC 열: 용량={c_cap}, soh={c_soh}, 시간={c_ts}, 조건={c_cond}, 누적충전량={c_q}")
for name, c in [("용량", c_cap), ("시간", c_ts), ("조건(CU 구분)", c_cond), ("누적 충전량", c_q)]:
    if c is None: notes.append(f"EOC '{name}' 열을 자동으로 못 찾음 -> K0_columns.txt 보고 확인 필요")

cond_values, rows = None, []
for f in eoc_files:
    cid, pid = cell_id(f)
    try: d = read_csv_any(f)
    except Exception as ex: notes.append(f"{cid} EOC 읽기 실패: {ex}"); continue
    if c_cond and cond_values is None: cond_values = d[c_cond].value_counts().head(10).to_dict()
    cu = d
    if c_cond:   # 조건 값에 'cu'가 들어 있는 행을 check-up으로 봄 (값 형식 미확인 -> 아래 보고에 값 목록 표시)
        m = d[c_cond].astype(str).str.lower().str.contains("cu")
        if m.any(): cu = d[m]
    cap = cu[c_cap].dropna() if c_cap else pd.Series(dtype=float)
    if len(cap) < 2: rows.append(dict(cell=cid, pid=pid, n_cu=len(cap))); continue
    soh = cap / cap.iloc[0]
    reach = bool((soh <= SOH_TH).any())
    t_h = (cu.loc[cap.index, c_ts] - cu.loc[cap.index[0], c_ts]) / 3600 if c_ts else pd.Series(np.nan, index=cap.index)
    t85 = float(t_h[soh <= SOH_TH].iloc[0]) if reach and c_ts else np.nan
    rows.append(dict(cell=cid, pid=pid, n_cu=len(cap), soh_last=float(soh.iloc[-1]), reach85=reach,
                     t85_h_approx=t85, has_time=c_ts is not None, has_q=c_q is not None,
                     soh_cap_last=float(cu[c_soh].dropna().iloc[-1]) if c_soh and cu[c_soh].notna().any() else np.nan))
E = pd.DataFrame(rows)
n_reach = int(E.get("reach85", pd.Series(dtype=bool)).fillna(False).sum())
print(f"[①] SOH {SOH_TH} 도달 셀: {n_reach} / {len(E)}  (판정 기준 10개 이상)")
print(f"    조건 열 값 상위: {cond_values}  <- 'CU'가 보이지 않으면 check-up 구분 방법 확인 필요")

# ---------- ③ EIS: 첫 check-up 상온 EIS ----------
eis_ok, eis_info = [], {}
if eis_files:
    e0 = read_csv_any(eis_files[0], nrows=5); ecols = list(e0.columns)
    c_ets = pick(ecols, "timestamp", avoid=("origin",)); c_f = pick(ecols, "freq")
    c_T = pick(ecols, "t_avg") or pick(ecols, "temp"); c_soc = pick(ecols, "soc")
    c_re = pick(ecols, "z_re_comp"); c_im = pick(ecols, "z_im_comp"); c_econd = pick(ecols, "cond")
    eis_info = dict(시간=c_ets, 주파수=c_f, 온도=c_T, SoC=c_soc, 실수부=c_re, 허수부=c_im, 조건=c_econd)
    print(f"[③] 자동으로 찾은 EIS 열: {eis_info}")
    for f in eis_files:
        cid, _ = cell_id(f)
        try: d = read_csv_any(f)
        except Exception as ex: notes.append(f"{cid} EIS 읽기 실패: {ex}"); continue
        if not (c_ets and c_T): eis_ok.append(dict(cell=cid, first_rt=np.nan)); continue
        first_day = d[c_ets].min() + 3 * 86400                        # 첫 check-up = 첫 측정 후 3일 이내 (가정)
        first = d[d[c_ets] <= first_day]
        rt = first[(first[c_T] >= RT_RANGE[0]) & (first[c_T] <= RT_RANGE[1])]
        eis_ok.append(dict(cell=cid, first_rt=len(rt) > 0,
                           socs=sorted(rt[c_soc].dropna().unique().tolist()) if c_soc and len(rt) else [],
                           fmin=float(d[c_f].min()) if c_f else np.nan, fmax=float(d[c_f].max()) if c_f else np.nan))
    for name, c in eis_info.items():
        if c is None: notes.append(f"EIS '{name}' 열을 자동으로 못 찾음 -> 확인 필요")
S = pd.DataFrame(eis_ok)
n_rt = int(S.first_rt.fillna(False).sum()) if len(S) else 0
print(f"[③] 첫 check-up에 상온 EIS가 있는 셀: {n_rt} / {len(S)}  (상온 = 평균 온도 {RT_RANGE[0]}~{RT_RANGE[1]} °C, 첫 check-up = 첫 측정 후 3일 이내 — 가정)")

# ---------- ④ 조건표 ----------
cond = E[["cell", "pid"]].copy()
if cfg_files:
    cfg = []
    for f in cfg_files:
        cid, _ = cell_id(f)
        try: d = read_csv_any(f); cfg.append(dict(cell=cid, **{c: d[c].iloc[0] for c in d.columns[:20]}))
        except Exception as ex: notes.append(f"{cid} CFG 읽기 실패: {ex}")
    if cfg: cond = cond.merge(pd.DataFrame(cfg), on="cell", how="left")
else:
    notes.append("CFG 파일이 없음 -> 조건(온도, C-rate, SoC 창)은 'Cycling Experiment Cell Overview' 엑셀에서 옮겨야 함")
cond = cond.merge(E.drop(columns=["pid"]), on="cell", how="left")
if len(S): cond = cond.merge(S, on="cell", how="left")
cond.to_csv(KOUT / "K0_conditions.csv", index=False, encoding="utf-8-sig")
print(f"[④] 조건표 저장: {KOUT / 'K0_conditions.csv'}")

# ---------- K0_audit.md 작성 ----------
verdict = "진행" if n_reach >= 10 else "보류 (정오 sync에서 기준선 조정 여부 논의)"
md = f"""# K0_audit — KIT 데이터 점검 결과 (트랙 B, 과업 B2)

> `stepK0_audit.py`가 자동으로 작성. **"확인 필요"는 사람이 확인해야 함.** 열 이름 전체 목록: `out/K/K0_columns.txt`

## 0. 내려받은 파일
- 폴더: `{KIT}`
- EOC {len(eoc_files)}개, EIS {len(eis_files)}개, CFG {len(cfg_files)}개, 남은 압축 파일 {len(archives)}개

## ① SOH {SOH_TH} 도달 셀 수
- SOH = 각 check-up 용량 / 첫 check-up 용량 (사전 등록과 동일). 용량 열: `{c_cap}`
- **도달 셀 수: {n_reach} / {len(E)}**  (기준 10개 이상)
- check-up 구분: 조건 열 `{c_cond}`의 값 중 'cu'가 들어간 행. 조건 열 값 상위: `{cond_values}`

## ② check-up별 시간·누적 충전량
- 시간 열: `{c_ts}` ({'있음' if c_ts else '없음 — 확인 필요'})
- 누적 충전량 열: `{c_q}` ({'있음 — EFC = 누적 Ah / (2 x 3 Ah)로 계산 가능' if c_q else '없음 — LOG_AGE 필요 여부 확인'})

## ③ 첫 check-up 상온 EIS
- **첫 check-up에 상온 EIS가 있는 셀: {n_rt} / {len(S)}**
- 자동으로 찾은 EIS 열: `{eis_info}`
- 가정: 상온 = 평균 온도 {RT_RANGE[0]}~{RT_RANGE[1]} °C, 첫 check-up = 첫 EIS 측정 후 3일 이내

## ④ 셀별 조건표
- `out/K/K0_conditions.csv` (셀, 파라미터 ID, cfg 열, check-up 수, 마지막 SOH, 0.85 도달, 첫 상온 EIS)

## 발견한 문제 / 확인 필요
""" + "\n".join(f"- {n}" for n in notes) + f"""
- 온도 이상 구간(2022-11-27~29, 2023-04 25 °C 풀) 표시: 아직 안 함 — 로더(B3)에서 처리

## 결론 (정오 sync용)
- KIT 사용: **{verdict}** — SOH {SOH_TH} 도달 {n_reach}셀
"""
(KOUT / "K0_audit.md").write_text(md, encoding="utf-8")
print(f"\n저장: {KOUT / 'K0_audit.md'}\n결론: KIT 사용 {verdict}")
if notes: print("확인 필요:", *notes, sep="\n  - ")
