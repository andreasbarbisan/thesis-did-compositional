clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global out "$root/checkpoints"

capture log close _all
log using "$out/logs/00e_testa_csdid_controles_municipais.log", text replace

use "$out/painel_completo.dta", clear
bysort year: egen double medprof_year_mean = mean(medprof)
bysort year: egen double medprof_year_sd = sd(medprof)
gen double medprof_std = (medprof - medprof_year_mean) / medprof_year_sd if !missing(medprof) & medprof_year_sd > 0
replace gvar = 0 if first_treat > 2018 & !missing(first_treat)
replace treated = 0 if first_treat > 2018 & !missing(first_treat)
keep if serie == "9EF" & co_comp == 1

local comp "pct_fem_imp pct_nee_imp idade_media_imp pct_preta_parda_imp pct_raca_declarada_imp pct_bolsa_fam_imp ln_alunos"
local mun "ln_pib_pc ln_pop prop_jovem_0_19"

di "=== DRDID original ==="
capture noisily csdid medprof_std medprof_base i.zona_id i.tipo_id `comp', ivar(codesc) time(year) gvar(gvar) method(dripw) notyet
di "rc original = " _rc

di "=== DRDID + municipais sem factor vars ==="
capture noisily csdid medprof_std medprof_base `mun' `comp', ivar(codesc) time(year) gvar(gvar) method(dripw) notyet
di "rc mun sem factor = " _rc

di "=== DRDID + municipais com factor vars ==="
capture noisily csdid medprof_std medprof_base i.zona_id i.tipo_id `mun' `comp', ivar(codesc) time(year) gvar(gvar) method(dripw) notyet
di "rc mun factor = " _rc

di "=== IPW + municipais ==="
capture noisily csdid medprof_std medprof_base `mun' `comp', ivar(codesc) time(year) gvar(gvar) method(ipw) notyet
di "rc ipw = " _rc

log close
