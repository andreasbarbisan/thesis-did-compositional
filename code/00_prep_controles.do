/*===========================================================================
  00_prep_controles.do - Prepara covariadas agregadas de matricula
  Andreas Azambuja Barbisan - TCC Insper 2026

  Papel atual no projeto:
    - Usar microdados de matricula como fonte de composicao agregada
    - Colapsar para escola-ano e escola-ano-serie
    - Nao tentar link individual aluno-nota

  Outputs:
    checkpoints/matricula_escola_ano.dta
    checkpoints/matricula_escola_ano_serie.dta

  Unidades:
    codesc x year
    codesc x year x serie, com serie em {9EF, EM3}
===========================================================================*/

clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global mat  "$root/Matrícula"
global out  "$root/checkpoints"

capture mkdir "$out"

di ""
di "=== PARTE 1: Matricula - covariadas agregadas ==="

tempfile emp_school emp_series
local first_school = 1
local first_series = 1

foreach ano in 2011 2012 2013 2014 2015 2016 2017 2018 {

    local arq "$mat/microdados_matricula_sp_`ano'_12.`ano'.csv"
    di "  `ano'...", _continue

    import delimited "`arq'", delimiter(";") encoding("UTF-16") ///
        clear varnames(1) case(lower)

    * Padronizar codigo da escola.
    capture confirm variable cd_escola
    if _rc {
        capture rename cod_escola cd_escola
        capture rename codesc cd_escola
    }

    * Filtrar rede estadual. PEI e politica estadual.
    keep if strpos(upper(nomedep), "ESTADUAL") > 0

    * Manter matriculas ativas quando o campo existir.
    capture confirm variable flag_sit_aluno
    if !_rc {
        capture destring flag_sit_aluno, replace force
        keep if missing(flag_sit_aluno) | flag_sit_aluno == 0
    }

    * Variaveis numericas basicas.
    foreach v in idade rendimento corraca end_zona serie turno durclasse flag_bolsa_fam {
        capture confirm variable `v'
        if !_rc {
            capture destring `v', replace force
        }
    }

    * Indicadores de composicao.
    gen byte um = 1
    gen byte fem = (upper(strtrim(sexo)) == "F") if !missing(sexo)
    gen byte masc = (upper(strtrim(sexo)) == "M") if !missing(sexo)

    gen byte nee = 0
    forvalues d = 1/10 {
        capture confirm variable def`d'
        if !_rc {
            capture destring def`d', replace force
            replace nee = 1 if def`d' > 0 & !missing(def`d')
        }
    }

    gen double idade_num = idade if !missing(idade)

    * CORRACA costuma seguir codificacao administrativa:
    * 1 branca, 2 preta, 3 parda, 4 amarela, 5 indigena, 6 nao declarada.
    gen byte raca_branca = (corraca == 1) if !missing(corraca)
    gen byte raca_preta  = (corraca == 2) if !missing(corraca)
    gen byte raca_parda  = (corraca == 3) if !missing(corraca)
    gen byte raca_pp     = inlist(corraca, 2, 3) if !missing(corraca)
    gen byte raca_decl   = inlist(corraca, 1, 2, 3, 4, 5) if !missing(corraca)

    gen byte bolsa_fam = (flag_bolsa_fam == 1) if !missing(flag_bolsa_fam)
    gen byte zona_rural = (end_zona == 1) if !missing(end_zona)

    gen byte turno_integral = .
    capture confirm variable ds_turma
    if !_rc {
        replace turno_integral = (strpos(upper(ds_turma), "INTEGRAL") > 0) if !missing(ds_turma)
    }
    replace turno_integral = 1 if turno == 6 & !missing(turno)
    replace turno_integral = 0 if missing(turno_integral) & !missing(turno)

    * Mapa simples para as series estimadas no SARESP agregado.
    gen str3 serie_saresp = ""
    replace serie_saresp = "9EF" if serie == 9
    replace serie_saresp = "EM3" if serie == 12

    gen year = `ano'

    keep cd_escola year serie_saresp um fem masc nee idade_num raca_branca raca_preta raca_parda ///
         raca_pp raca_decl bolsa_fam zona_rural turno_integral
    drop if missing(cd_escola)

    tempfile raw_`ano'
    save `raw_`ano'', replace

    * Escola-ano: diagnostico geral da escola.
    preserve
        collapse ///
            (count) n_alunos=um ///
            (mean) pct_fem=fem pct_masc=masc pct_nee=nee idade_media=idade_num ///
                   pct_branca=raca_branca pct_preta=raca_preta pct_parda=raca_parda ///
                   pct_preta_parda=raca_pp pct_raca_declarada=raca_decl ///
                   pct_bolsa_fam=bolsa_fam pct_rural=zona_rural pct_integral=turno_integral, ///
            by(cd_escola year)
        rename cd_escola codesc

        if `first_school' {
            save `emp_school', replace
            local first_school = 0
        }
        else {
            append using `emp_school'
            save `emp_school', replace
        }
    restore

    * Escola-ano-serie: covariadas principais da estimacao.
    keep if inlist(serie_saresp, "9EF", "EM3")
    collapse ///
        (count) n_alunos=um ///
        (mean) pct_fem=fem pct_masc=masc pct_nee=nee idade_media=idade_num ///
               pct_branca=raca_branca pct_preta=raca_preta pct_parda=raca_parda ///
               pct_preta_parda=raca_pp pct_raca_declarada=raca_decl ///
               pct_bolsa_fam=bolsa_fam pct_rural=zona_rural pct_integral=turno_integral, ///
        by(cd_escola year serie_saresp)
    rename cd_escola codesc
    rename serie_saresp serie

    if `first_series' {
        save `emp_series', replace
        local first_series = 0
    }
    else {
        append using `emp_series'
        save `emp_series', replace
    }

    di " OK"
}

use `emp_school', clear
label var codesc             "Codigo da escola"
label var year               "Ano"
label var n_alunos           "Matriculas estaduais ativas"
label var pct_fem            "Fracao de alunas"
label var pct_masc           "Fracao de alunos homens"
label var pct_nee            "Fracao com necessidade especial registrada"
label var idade_media        "Idade media dos matriculados"
label var pct_branca         "Fracao branca"
label var pct_preta          "Fracao preta"
label var pct_parda          "Fracao parda"
label var pct_preta_parda    "Fracao preta ou parda"
label var pct_raca_declarada "Fracao com raca/cor declarada"
label var pct_bolsa_fam      "Fracao com Bolsa Familia"
label var pct_rural          "Fracao com endereco rural"
label var pct_integral       "Fracao em turno integral"
sort codesc year
save "$out/matricula_escola_ano.dta", replace

use `emp_series', clear
label var codesc             "Codigo da escola"
label var year               "Ano"
label var serie              "Serie SARESP aproximada pela matricula"
label var n_alunos           "Matriculas estaduais ativas na serie"
label var pct_fem            "Fracao de alunas na serie"
label var pct_masc           "Fracao de alunos homens na serie"
label var pct_nee            "Fracao com necessidade especial registrada na serie"
label var idade_media        "Idade media dos matriculados na serie"
label var pct_branca         "Fracao branca na serie"
label var pct_preta          "Fracao preta na serie"
label var pct_parda          "Fracao parda na serie"
label var pct_preta_parda    "Fracao preta ou parda na serie"
label var pct_raca_declarada "Fracao com raca/cor declarada na serie"
label var pct_bolsa_fam      "Fracao com Bolsa Familia na serie"
label var pct_rural          "Fracao com endereco rural na serie"
label var pct_integral       "Fracao em turno integral na serie"
sort codesc year serie
save "$out/matricula_escola_ano_serie.dta", replace

di ""
di "Salvo: checkpoints/matricula_escola_ano.dta"
di "Salvo: checkpoints/matricula_escola_ano_serie.dta"
di ""
di "Covariadas por serie geradas:"
describe

di ""
di "=== Concluido. Proximo passo: 01_prep.do e depois 01b_merge_controles.do ==="
