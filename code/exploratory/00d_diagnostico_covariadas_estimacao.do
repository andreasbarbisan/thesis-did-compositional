clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global out "$root/checkpoints"

capture log close _all
log using "$out/logs/00d_diagnostico_covariadas_estimacao.log", text replace

use "$out/painel_completo.dta", clear

bysort year: egen double medprof_year_mean = mean(medprof)
bysort year: egen double medprof_year_sd = sd(medprof)
gen double medprof_std = (medprof - medprof_year_mean) / medprof_year_sd if !missing(medprof) & medprof_year_sd > 0
replace gvar = 0 if first_treat > 2018 & !missing(first_treat)
replace treated = 0 if first_treat > 2018 & !missing(first_treat)

keep if serie == "9EF" & co_comp == 1

foreach v in medprof_std medprof_base zona_id tipo_id ln_pib_pc ln_pop prop_jovem_0_19 pct_fem_imp pct_nee_imp idade_media_imp pct_preta_parda_imp pct_raca_declarada_imp pct_bolsa_fam_imp ln_alunos {
    qui count if missing(`v')
    di "`v' missing: " r(N) " / " _N
    tab treated if missing(`v')
}

gen byte complete_old = !missing(medprof_std, medprof_base, zona_id, tipo_id, pct_fem_imp, pct_nee_imp, idade_media_imp, pct_preta_parda_imp, pct_raca_declarada_imp, pct_bolsa_fam_imp, ln_alunos)
gen byte complete_new = !missing(medprof_std, medprof_base, zona_id, tipo_id, ln_pib_pc, ln_pop, prop_jovem_0_19, pct_fem_imp, pct_nee_imp, idade_media_imp, pct_preta_parda_imp, pct_raca_declarada_imp, pct_bolsa_fam_imp, ln_alunos)

tab complete_old treated
tab complete_new treated
tab gvar if complete_new

log close
