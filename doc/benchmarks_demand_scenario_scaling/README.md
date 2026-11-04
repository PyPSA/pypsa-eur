# Benchmark TSSB: tempo e memoria al variare degli scenari

Questa procedura documenta i grafici generati dai risultati già presenti nel repository, senza rilanciare le ottimizzazioni. I comandi seguenti vanno eseguiti dalla radice del repository, con l'ambiente Python/Snakemake del progetto attivo. Per la raccolta e i grafici è stato usato `/home/pampado/miniconda3/envs/pypsa-eur/bin/python`.

## Casi confrontati

Prefisso comune: `results/capacity-expansion-tssb-power-IT-noUC/`.

| Cartella | Numero di scenari di domanda |
| --- | ---: |
| `th_24c24_demand_20` | 20 |
| `th_24c24_demand_50` | 50 |
| `th_24c24_demand_100` | 100 |
| `th_24c24_demand_200` | 200 |

Si confrontano separatamente reti a **5 nodi** e a **20 nodi** (`clusters`), mantenendo la risoluzione temporale `24c24` (24 periodi di 24 ore). L'asse orizzontale rappresenta il numero di scenari stocastici; i nodi rimangono costanti in ciascuna figura.

| Identificativo nei file | Etichetta nei grafici (`solver_options`) |
| --- | --- |
| `gurobi` | `gurobi-default` |
| `smspp` | `smspp-tssb` |
| `smspp_2` | `smspp-tssb-bds` |
| `smspp_3` | `smspp-tssb-ip-par` |

La raccolta aggiornata il **4 novembre 2026** contiene 32 righe: 4 numerosità × 2 dimensioni di rete × 4 identificativi di solver. Su indicazione esplicita dell’utente, `smspp_3` è identificato come `smspp-tssb-ip-par` in tutti i casi, inclusa la run fallita a 20 nodi/200 scenari. Questa assegnazione sostituisce le etichette nei metadati delle reti e nei config di ripiego.

## Rule e configurazioni per generare le esecuzioni

Il target aggregatore corretto è **`compare_solver_elec_tssb_outputs`**, definito in [rules/collect.smk](../../rules/collect.smk). La catena delle dipendenze comprende:

1. `build_stochastic_elec_network`, in [rules/solve_stochastic.smk](../../rules/solve_stochastic.smk), che costruisce la rete stocastica a partire dalla rete elettrica preparata.
2. `solve_network`, in [rules/solve_electricity.smk](../../rules/solve_electricity.smk), che risolve la rete con ciascun solver configurato, tramite `scripts/solve_network_supervisor.py`.
3. `make_solver_comparison_elec_tssb`, in [rules/postprocess.smk](../../rules/postprocess.smk), che legge reti risolte e benchmark ed esegue `scripts/compare_solver_results.py`.
4. `compare_solver_elec_tssb_outputs`, che richiede i CSV di confronto per i casi e i cluster selezionati.

La configurazione di istanza è [config/instances/tssb/power-IT-noUC.yaml](../../config/instances/tssb/power-IT-noUC.yaml), integrata dai default di `config/config.default.yaml` e dagli override per caso. I campi rilevanti sono:

```yaml
run:
  prefix: capacity-expansion-tssb-power-IT-noUC
  name: th_24c24_demand_20
  scenarios:
    enable: true
    file: config/instances/scenarios/tssb/scenarios-power.yaml
scenario:
  clusters: [5, 20]
  opts: ['']
countries: [IT]
stochastic_scenarios:
  enable: true
  file: config/instances/scenarios/tssb/stochastic-demand_20.yaml
clustering:
  temporal:
    resolution_elec: 24c24
solving:
  solver:
    name: [gurobi, smspp, smspp, smspp]
    options: [gurobi-default, smspp-tssb, smspp-tssb-bds, smspp-tssb-ip-par]
```

Questo estratto illustra i parametri necessari; non sostituisce il file completo. Nelle configurazioni salvate sono presenti anche `skip_iterations: true`, `transmission_losses: 0` e `linearized_unit_commitment: false`.

Occorre distinguere i **casi del workflow** (`run.scenarios.file`, che assegna override a nomi come `th_24c24_demand_20`) dagli **scenari stocastici interni alla rete** (`stochastic_scenarios.file`). Questi ultimi sono definiti nei file:

```text
config/instances/scenarios/tssb/stochastic-demand_20.yaml
config/instances/scenarios/tssb/stochastic-demand_50.yaml
config/instances/scenarios/tssb/stochastic-demand_100.yaml
config/instances/scenarios/tssb/stochastic-demand_200.yaml
```

Gli scenari applicano una modifica a `Load.p_set` tramite `modify_components`, con operazione `scale`, distribuzione normale, media 1, deviazione standard 0.1, estremi 0.8 e 1.2 e seed 123. Le probabilità sono uniformi (`1/N`).

**Stato dei config al 4 novembre 2026:** il file di istanza attuale seleziona `th_24c24_demand_20` e il solo cluster 20; inoltre `scenarios-power.yaml` contiene il caso `th_24c24_demand_20`, ma non i casi `th_24c24_demand_50`, `_100` e `_200`. Non basta quindi lanciare il config attuale per riprodurre tutte e quattro le cartelle. Per ricrearle, aggiungere al file degli override i tre casi mancanti, copiando il blocco `th_24c24_demand_20`, aggiornando il nome e il percorso `stochastic-demand_N.yaml`, e impostando `parameters: [demand]` per le opzioni SMS++ utilizzate.

Poi selezionare nel config di istanza i quattro nomi in `run.name` e `scenario.clusters: [5, 20]`. Il comando per richiedere l'intera catena è:

```bash
snakemake --cores 32 \
  --configfile config/instances/tssb/power-IT-noUC.yaml \
  compare_solver_elec_tssb_outputs
```

Questo è un comando di riproduzione, non un comando rilanciato per realizzare i grafici. I parametri effettivamente esportati dalle esecuzioni storiche sono nei file `configs/config.base_s_<C>_elec__<SOLVER>.yaml` di ciascuna cartella; vanno consultati per confronti precisi, insieme ai file SMS++ referenziati. Per esempio, il config salvato a 5 nodi/20 scenari indica `gurobi-default.threads: 32`, `TimeLimit: 43200`, `smspp-tssb.configfile: config/smspp/TSSBlock/TSSBCfg_grb.txt` e `smspp-tssb-bds.configfile: config/smspp/TSSBlock/TSSBSCfg-BDS.txt`. I file di configurazione esterni possono essere cambiati dopo la run.

## File prodotti dalle esecuzioni

Indicando con `N` il numero di scenari, `C` i nodi e `SOLVER` uno fra `gurobi`, `smspp`, `smspp_2`, nella cartella `results/capacity-expansion-tssb-power-IT-noUC/th_24c24_demand_N/` si trovano:

```text
networks/base_s_C_elec__SOLVER.nc
configs/config.base_s_C_elec__SOLVER.yaml
benchmarks/solve_network/base_s_C_elec__SOLVER
benchmarks/solve_network/base_s_C_elec__SOLVER.smspp
benchmarks/solve_network/base_s_C_elec__SOLVER.failed
logs/solve_network/base_s_C_elec__SOLVER_python.log
csvs/solver_comparison/summary_s_C_elec__tssb.csv
csvs/solver_comparison/benchmarks_s_C_elec__tssb.csv
```

I suffissi `.smspp` riguardano i tempi interni SMS++; i file `.failed` non vuoti contengono le misure di esecuzioni fallite. Snakemake può creare placeholder `.failed` vuoti anche per run riuscite: il collector distingue i due casi.

I grafici qui documentati leggono i benchmark di **`solve_network`**, attraverso il collector seguente. Non usano come input i CSV `csvs/solver_comparison/benchmarks_*.csv`, né i benchmark della costruzione della rete o del postprocessing. Il tempo rappresenta quindi la fase della rule di solve, non l'intero workflow.

## Raccolta dei dati per i grafici

Script: [scripts/collect_solve_benchmarks.py](../../scripts/collect_solve_benchmarks.py).

```bash
python scripts/collect_solve_benchmarks.py \
  --folder \
    capacity-expansion-tssb-power-IT-noUC/th_24c24_demand_20 \
    capacity-expansion-tssb-power-IT-noUC/th_24c24_demand_50 \
    capacity-expansion-tssb-power-IT-noUC/th_24c24_demand_100 \
    capacity-expansion-tssb-power-IT-noUC/th_24c24_demand_200 \
  --output results/capacity-expansion-tssb-power-IT-noUC/benchmarks_demand_scenario_scaling/solve_benchmarks.csv
```

Il collector legge l'ultima riga di ciascun benchmark TSV, aggiunge i metadati delle reti NetCDF e gli esiti dai log Python. Per le run fallite usa il benchmark `.failed`; se manca il config specifico può recuperare l'identità del solver da un config della stessa variante presente sotto il prefisso comune. Questi metadati di ripiego non dimostrano che tutti i parametri delle due esecuzioni coincidano.

| Campo originale / origine | Campo usato nel plot | Significato |
| --- | --- | --- |
| `s`, benchmark principale | `runtime_s` | Tempo trascorso in secondi della rule di solve; linea continua e punti. |
| `max_rss`, benchmark principale | `memory_gib = max_rss / 1024` | Picco della memoria residente registrato dal benchmark, da MiB a GiB. |
| `s`, file `.smspp` | `smspp_optimization_s` | Durata del sottoprocesso SMS++; linea puntinata. |
| `computational_time`, file `.smspp` | `smspp_computational_s` | Tempo computazionale riportato da SMS++; linea tratteggiata. |
| Dimensione `scenario` della rete NetCDF | `scenarios` | Numero di scenari; se mancante, dedotto dal suffisso numerico di `case`. |
| Nome `base_s_C_elec__SOLVER` | `clusters` | Numero di nodi usato per filtrare ciascun grafico. |
| Metadati config/rete | `solver_options` | Etichetta e colore della curva del solver. |
| Percorso della cartella | `case`, poi `periods`, `hours_per_period`, `case_variant` | Risoluzione temporale e variante (`demand`). |
| Log Python / suffisso `.failed` | `solve_status`, `termination_condition`, poi `outcome` | Stile dei punti in funzione dell'esito. |
| Config e log | `time_limit_s` | Limite usato come ordinata quando il testo dell'esito segnala esplicitamente un time limit. |

`cpu_time`, `h:m:s`, `max_vms`, `max_uss`, `max_pss`, I/O e carico medio sono conservati nel CSV raccolto, ma non determinano le ordinate di questi grafici. `objective` e `objective_constant` vengono letti dalla rete e preparati dal codice comune, ma qui non sono rappresentati.

Dopo la raccolta, applicare la correzione dell’identità di `smspp_3` confermata dall’utente, prima di generare i grafici:

```bash
python - <<'PY'
from pathlib import Path
import pandas as pd
path = Path("results/capacity-expansion-tssb-power-IT-noUC/benchmarks_demand_scenario_scaling/solve_benchmarks.csv")
data = pd.read_csv(path)
data.loc[data["solver"].eq("smspp_3"), "solver_options"] = "smspp-tssb-ip-par"
data.to_csv(path, index=False)
PY
```

## Generazione dei grafici

Script: [scripts/plot_benchmark_scenarios.py](../../scripts/plot_benchmark_scenarios.py), che riusa preparazione dati, colori, linee e legenda di [scripts/plot_solve_benchmarks.py](../../scripts/plot_solve_benchmarks.py).

```bash
python scripts/plot_benchmark_scenarios.py \
  results/capacity-expansion-tssb-power-IT-noUC/benchmarks_demand_scenario_scaling/solve_benchmarks.csv \
  --nodes 5

python scripts/plot_benchmark_scenarios.py \
  results/capacity-expansion-tssb-power-IT-noUC/benchmarks_demand_scenario_scaling/solve_benchmarks.csv \
  --nodes 20
```

La preparazione comune elimina righe con nodi o runtime non validi e aggrega eventuali ripetizioni tramite mediana, separando variante, numero di scenari, solver, nodi, risoluzione ed esito. Poi si filtra il numero di nodi richiesto. Lo script richiede un'unica variante e un'unica risoluzione temporale; ordina i punti per numero di scenari. Il tempo usa l'asse verticale logaritmico; la memoria usa l'asse verticale lineare e non mostra i tempi interni SMS++.

Gli output sono in:

```text
results/capacity-expansion-tssb-power-IT-noUC/benchmarks_demand_scenario_scaling/
├── solve_benchmarks.csv
└── plots/
    ├── benchmark_scenario_scaling_5_nodes.png
    ├── benchmark_scenario_scaling_5_nodes.pdf
    ├── benchmark_memory_scenario_scaling_5_nodes.png
    ├── benchmark_memory_scenario_scaling_5_nodes.pdf
    ├── benchmark_scenario_scaling_5_nodes.csv
    ├── benchmark_scenario_scaling_20_nodes.png
    ├── benchmark_scenario_scaling_20_nodes.pdf
    ├── benchmark_memory_scenario_scaling_20_nodes.png
    ├── benchmark_memory_scenario_scaling_20_nodes.pdf
    └── benchmark_scenario_scaling_20_nodes.csv
```

I CSV sotto `plots/` contengono i dati preparati e filtrati, inclusi `runtime_s`, `memory_gib`, tempi SMS++, `scenarios`, `solver_options` e `outcome`.

## Lettura degli esiti e limite delle figure attuali

Per 5 nodi tutte le 16 righe sono classificate dal codice come `solved`: ciò indica che non sono stati rilevati i marcatori di fallimento/time limit/infeasibilità, non certifica da solo l'ottimalità di ogni soluzione.

Per **20 nodi e 200 scenari**, le tre esecuzioni SMS++ provengono da file `.failed` e riportano `solve_status=error`, `termination_condition=exit_code_1`. Le misure aggiornate sono:

| Identificativo | Etichetta raccolta | Tempo [s] | Picco RSS [GiB] |
| --- | --- | ---: | ---: |
| `smspp` | `smspp-tssb` | 4470.67 | 152.40 |
| `smspp_2` | `smspp-tssb-bds` | 4557.39 | 85.99 |
| `smspp_3` | `smspp-tssb-ip-par` (indicazione utente) | 1866.23 | 62.28 |

Non sono tempi di soluzione completata.

Il codice comune attuale assegna `outcome=time_limit` anche a ogni benchmark `.failed`, indipendentemente dalla causa. Per questo nelle figure i tre punti appaiono con una X e legenda **Time limit**, senza collegamento alla curva. **L'etichetta non prova che sia stato raggiunto un limite temporale:** in questi dati l'esito è un errore `exit_code_1`. Per stabilirne la causa occorre analizzare i log della singola run. I grafici e i CSV sono stati rigenerati il 4 novembre 2026 mantenendo questo comportamento del codice comune.
