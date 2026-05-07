/*===========================================================================
  00c_inspeciona_controles_municipais.do
  Inspeciona a base municipal externa para validar chaves e nomes de variaveis.
===========================================================================*/

clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global out  "$root/checkpoints"
global ctrl "C:/Users/andreas/OneDrive - Insper - Instituto de Ensino e Pesquisa/Insper/Pesquisa/IC/Data/Controles/ControlesTCCAdriano.dta"

capture log close _all
log using "$out/logs/00c_inspeciona_controles_municipais.log", text replace

use "$ctrl", clear

describe
codebook, compact

ds, has(type numeric)
local numvars `r(varlist)'

di ""
di "=== Variaveis numericas candidatas ==="
di "`numvars'"

foreach pat in cod ibge mun ano year pib pop jovem idade {
    di ""
    di "=== Busca por padrao: `pat' ==="
    capture ds *`pat'*
    if !_rc di "`r(varlist)'"
}

log close
