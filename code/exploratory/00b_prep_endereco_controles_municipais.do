/*===========================================================================
  00b_prep_endereco_controles_municipais.do
  Prepara cadastro de escolas e controles municipais para o painel PEI.

  Inputs:
    TCC v2 Micro 4/bruto/ENDERECO_ESCOLAS_062025.csv
    ControlesTCCAdriano.dta

  Outputs:
    checkpoints/endereco_escolas.dta
    checkpoints/controles_municipais.dta
    checkpoints/diagnostico_merge_municipal.dta/.csv
===========================================================================*/

clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global bruto "$root/TCC v2 Micro 4/bruto"
global out  "$root/checkpoints"
global ctrl "C:/Users/andreas/OneDrive - Insper - Instituto de Ensino e Pesquisa/Insper/Pesquisa/IC/Data/Controles/ControlesTCCAdriano.dta"

capture mkdir "$out"
capture mkdir "$out/logs"
capture log close _all
log using "$out/logs/00b_prep_endereco_controles_municipais.log", text replace

*===========================================================================
* 1. Cadastro de escolas
*===========================================================================

import delimited "$bruto/ENDERECO_ESCOLAS_062025.csv", ///
    delimiter(";") encoding("utf-8") clear case(lower) varnames(1)

keep nomedep de mun cd_ibge distr cod_esc nomesc codsit tipoesc zona ///
    ds_latitude ds_longitude

rename cod_esc codesc
rename cd_ibge cd_ibge7
rename ds_latitude latitude
rename ds_longitude longitude

foreach v in latitude longitude {
    replace `v' = subinstr(`v', ",", ".", .)
    destring `v', replace force
}

gen long cd_ibge6 = floor(cd_ibge7 / 10)

foreach v in de mun distr nomesc {
    replace `v' = strtrim(upper(`v'))
}

duplicates drop codesc, force

label var codesc    "Codigo da escola"
label var cd_ibge7  "Codigo IBGE municipal com 7 digitos"
label var cd_ibge6  "Codigo IBGE municipal com 6 digitos"
label var de        "Diretoria de ensino"
label var mun       "Municipio"
label var distr     "Distrito"
label var latitude  "Latitude da escola"
label var longitude "Longitude da escola"

save "$out/endereco_escolas.dta", replace

*===========================================================================
* 2. Controles municipais
*===========================================================================

use "$ctrl", clear

keep if inrange(year, 2011, 2018)
keep ibge year pib pop total_0_4 total_5_9 total_10_14 total_15_19

rename ibge cd_ibge6
gen double pib_pc = pib / pop if pop > 0
gen double ln_pib_pc = log(pib_pc) if pib_pc > 0
gen double ln_pop = log(pop) if pop > 0
gen double prop_jovem_0_19 = total_0_4 + total_5_9 + total_10_14 + total_15_19

label var pib              "PIB municipal em R$ correntes"
label var pop              "Populacao municipal"
label var pib_pc           "PIB municipal per capita"
label var ln_pib_pc        "Log do PIB per capita municipal"
label var ln_pop           "Log da populacao municipal"
label var prop_jovem_0_19  "Proporcao da populacao de 0 a 19 anos"

keep cd_ibge6 year pib pop pib_pc ln_pib_pc ln_pop prop_jovem_0_19
duplicates drop cd_ibge6 year, force

save "$out/controles_municipais.dta", replace

*===========================================================================
* 3. Diagnostico do merge potencial
*===========================================================================

use "$out/endereco_escolas.dta", clear
expand 8
bysort codesc: gen year = 2010 + _n
keep if inrange(year, 2011, 2018)
merge m:1 cd_ibge6 year using "$out/controles_municipais.dta", keep(master match) nogen
gen byte tem_controle_mun = !missing(ln_pib_pc, ln_pop, prop_jovem_0_19)
collapse (count) n_escola_ano=codesc (sum) n_com_controle=tem_controle_mun, by(year)
gen cobertura = n_com_controle / n_escola_ano

save "$out/diagnostico_merge_municipal.dta", replace
export delimited using "$out/diagnostico_merge_municipal.csv", replace

di ""
di "Preparo municipal concluido:"
di "  - $out/endereco_escolas.dta"
di "  - $out/controles_municipais.dta"
di "  - $out/diagnostico_merge_municipal.csv"

log close
