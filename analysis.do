
version 18.0
clear all
set more off
set maxvar 12000
set scheme s1color

global root   "/Users/micoria/Documents/research_projects/cognitive/B960-CLHLS数据清洗/CLHLS"
global work   "$root/Working_data"
global tables "$root/Tables"
global figs   "$root/Figures"
cap mkdir "$tables"
cap mkdir "$figs"

cap which esttab
if _rc ssc install estout, replace
cap which coefplot
if _rc ssc install coefplot, replace


* --- 生化表：去重 → 固化 ---
use "$work/biomarker_2014_clean.dta", clear
destring id, replace force
format id %20.0f

quietly duplicates report id
if r(N)>0 & r(unique_value)<r(N) {
    di as result "生化表存在重复ID，聚合（连续=mean，四分类=firstnm）"
    preserve
    collapse (mean) alb glu bun crea cho tg gsp crphs ua hdlc sod mda vd3 vitb12 ///
                     ualb ucr ualbbyucr wbc lymph lymph_a rbc hgb hct mcv mch mchc ///
                     plt mpv pdw pct non_hdl aip bio_abn_burden ///
                     (firstnm) alb_cat glu_cat bun_cat crea_cat cho_cat tg_cat gsp_cat ua_cat ///
                               hdlc_cat sod_cat wbc_cat lymph_cat lymph_a_cat rbc_cat hgb_cat hct_cat mcv_cat mch_cat mchc_cat ///
                               plt_cat mpv_cat pdw_cat pct_cat, by(id)
    save "$work/__biomarker_ready_2014.dta", replace
    restore
}
else {
    save "$work/__biomarker_ready_2014.dta", replace
}



use "$work/clhls08_18.dta", clear
destring id, replace force
format id %20.0f

gen in2014 = !missing(mmse_14) | !missing(trueage_14) | !missing(adl_d_14)
keep if in2014==1
drop in2014

quietly duplicates report id
if r(N)>0 & r(unique_value)<r(N) {
    di as result "2014 子样本存在重复ID，保留第一条"
    bys id: keep if _n==1
}
isid id, sort

* --- 合并 ---
merge 1:1 id using "$work/__biomarker_ready_2014.dta"

di as text "—— 合并结果 ——"
tab _merge
count if _merge==3
local matched = r(N)
count
local total2014 = r(N)
di as result "2014 子样本: `total2014'；匹配生化: `matched'；匹配率=" %6.2f (100*`matched'/`total2014') "%"

drop if _merge==2
drop _merge
save "$work/biomarker_merged_2014.dta", replace

*以14年为基线
gen age14 = trueage_14
replace age14 = trueage_11 if missing(age14) & !missing(trueage_11)
replace age14 = trueage_08 if missing(age14) & !missing(trueage_08)
label var age14 "年龄(2014为主)"

gen edu14 = edu_14
replace edu14 = edu_08 if missing(edu14) & !missing(edu_08)
label var edu14 "受教育年限(2014为主)"

gen smoke14 = smoke_14
replace smoke14 = smoke_11 if missing(smoke14)
replace smoke14 = smoke_08 if missing(smoke14)
label var smoke14 "是否吸烟(2014优先)"

gen drink14 = drink_14
replace drink14 = drink_11 if missing(drink14)
replace drink14 = drink_08 if missing(drink14)
label var drink14 "是否饮酒(2014优先)"

gen exercise14 = exercise_14
replace exercise14 = exercise_11 if missing(exercise14)
replace exercise14 = exercise_08 if missing(exercise14)
label var exercise14 "是否锻炼(2014优先)"

gen adl14 = adl_d_14
label var adl14 "ADL残疾(2014)"


foreach v in alb_cat glu_cat bun_cat crea_cat cho_cat tg_cat gsp_cat ua_cat hdlc_cat ///
             sod_cat wbc_cat lymph_cat lymph_a_cat rbc_cat hgb_cat hct_cat mcv_cat  ///
             mch_cat mchc_cat plt_cat mpv_cat pdw_cat pct_cat {
    cap fvset base 3 `v'
}

local catvars alb_cat glu_cat bun_cat crea_cat cho_cat tg_cat gsp_cat ua_cat hdlc_cat ///
              sod_cat wbc_cat lymph_cat lymph_a_cat rbc_cat hgb_cat hct_cat mcv_cat  ///
              mch_cat mchc_cat plt_cat mpv_cat pdw_cat pct_cat
foreach v of local catvars {
    cap drop abn_`v'
    gen abn_`v' = inlist(`v',1,2) if !missing(`v')
}
cap drop bio_abn_burden
egen bio_abn_burden = rowtotal(abn_*), missing
label var bio_abn_burden "生化异常负荷(异常项总数)"

foreach v in glu tg hdlc ua cho {
    cap drop `v'_q
    xtile `v'_q = `v', nq(4)
    label define q4 1 "Q1最低" 2 "Q2" 3 "Q3" 4 "Q4最高"
    label values `v'_q q4
}


egen mmse_count = rownonmiss(mmse_08 mmse_11 mmse_14 mmse_18)
gen has_bio = !missing(glu) | !missing(tg) | !missing(hdlc) | !missing(cho) | !missing(ua)


* 主分析 1：纵向认知（LME：随机截距+斜率）
preserve
keep if !missing(mmse_14) & mmse_count>=2 & has_bio
keep id sex age14 edu14 smoke14 drink14 exercise14 adl14 ///
     mmse_08 mmse_11 mmse_14 mmse_18 ///
     glu_cat hdlc_cat tg_cat cho_cat ua_cat bio_abn_burden

capture confirm variable mmse_08
if !_rc rename mmse_08 mmse_8

reshape long mmse_, i(id) j(year)
rename mmse_ mmse

gen time_rel14 = .
replace time_rel14 = -6 if year==8
replace time_rel14 = -3 if year==11
replace time_rel14 =  0 if year==14
replace time_rel14 =  4 if year==18
label var time_rel14 "相对2014的年份"

mixed mmse i.glu_cat##c.time_rel14 i.hdlc_cat##c.time_rel14 i.tg_cat##c.time_rel14 ///
      c.bio_abn_burden##c.time_rel14 ///
      c.age14 i.sex c.edu14 i.smoke14 i.drink14 i.exercise14 i.adl14 ///
      || id: time_rel14, covariance(unstructured) vce(robust)
estimates store LME1
esttab LME1 using "$tables/LME1_mmse_decline_2014centered.rtf", ///
      se star(* 0.05 ** 0.01 *** 0.001) label replace

* 轨迹图（按 TG 分类）
tempfile tmp
save `tmp'
collapse (mean) mmse , by(time_rel14 tg_cat)
twoway (line mmse time_rel14 if tg_cat==1, lp(solid) ms(o)) ///
       (line mmse time_rel14 if tg_cat==2, lp(dash)  ms(o)) ///
       (line mmse time_rel14 if tg_cat==3, lp(dot)   ms(o)), ///
       legend(order(1 "TG>上限" 2 "TG<下限" 3 "TG正常") pos(6) ring(0)) ///
       xtitle("相对2014的年份") ytitle("MMSE均值") ///
       title("MMSE轨迹（按TG分类）")
graph export "$figs/mmse_traj_by_TG.png", replace
use `tmp', clear
restore


* 主分析 2：2014–2018 死亡（KM + Cox；中点近似）
preserve
keep if dth11_14==0   // 进入2014存活

* 无具体死亡时间：死亡=2年中点；其余(存活/失访)=4年右删
gen followup_years = cond(dth14_18==1, 2, 4)
label var followup_years "随访时长(年)：死亡=2；其余=4（近似）"
gen died_1418 = (dth14_18==1)
label var died_1418 "2014–2018期间死亡(1=是)"

stset followup_years, failure(died_1418==1)

* Cox 回归
stcox i.glu_cat i.hdlc_cat i.tg_cat i.cho_cat i.ua_cat ///
      c.bio_abn_burden c.age14 i.sex c.edu14 i.smoke14 i.drink14 i.exercise14 i.adl14, vce(robust)
estimates store COX1
esttab COX1 using "$tables/COX1_mortality_2014to2018.rtf", eform b(%6.2f) ci(%6.2f) ///
      star(* 0.05 ** 0.01 *** 0.001) label replace

* KM 生存曲线（异常负荷三分位）
xtile burden_ter = bio_abn_burden if !missing(bio_abn_burden), nq(3)
label define tert 1 "低" 2 "中" 3 "高", replace
label values burden_ter tert
sts graph, by(burden_ter) risktable tmax(4) ///
    title("2014生化异常负荷与生存（至2018；死亡时间=2年近似）") ///
    xtitle("年") ytitle("生存概率")
graph export "$figs/survival_by_burden_tertiles.png", replace
restore

* 主分析 3：2018 认知障碍（Logit）
preserve
keep if !missing(mmse_d_18)
logit mmse_d_18 i.glu_cat i.hdlc_cat i.tg_cat i.cho_cat i.ua_cat ///
     c.bio_abn_burden c.age14 i.sex c.edu14 i.smoke14 i.drink14 i.exercise14 i.adl14, vce(robust)
estimates store LOGIT1
esttab LOGIT1 using "$tables/LOGIT1_mmseD18.rtf", eform b(%6.2f) ci(%6.2f) ///
      star(* 0.05 ** 0.01 *** 0.001) label replace

coefplot LOGIT1, keep(*glu_cat *hdlc_cat *tg_cat *cho_cat *ua_cat bio_abn_burden) ///
    eform xline(1) ciopts(recast(rcap)) vertical ///
    title("2014生化与2018认知障碍（OR）") ///
    ylabel(, angle(0)) xlabel(, angle(45))
graph export "$figs/forest_OR_biomarkers_mmseD18.png", replace
restore

* 敏感性分析
 *1) 生存：剔除2018失访者再跑（可能 ↑ 事件比例）
preserve
keep if dth11_14==0 & inlist(dth14_18,0,1)
gen followup_years = cond(dth14_18==1, 2, 4)
gen died_1418 = (dth14_18==1)
stset followup_years, failure(died_1418)
stcox i.glu_cat i.hdlc_cat i.tg_cat i.cho_cat i.ua_cat ///
      c.bio_abn_burden c.age14 i.sex c.edu14 i.smoke14 i.drink14 i.exercise14 i.adl14, vce(robust)
restore


di as text "== 全部完成。表：$tables ；图：$figs。=="

