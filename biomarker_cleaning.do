
clear all
set more off

cd "/Users/micoria/Documents/research_projects/cognitive/B960-CLHLS数据清洗/CLHLS/Raw_data/clhls/"
import delimited "biomarker_dataset_CLHLS_2014-1.tab", delimiter(tab) varn(1) stringcols(_all) case(preserve) clear

capture confirm variable id
if _rc==0 {
    replace id = subinstr(id, " ", "", .)
    destring id, replace force
    order id, first
}

* 将所有生化变量转为数值，负值与明显错误设为缺失
local biomvars alb glu bun crea cho tg gsp crphs ua hdlc sod mda vd3 vitb12 ///
               ualb ucr ualbbyucr wbc lymph lymph_a rbc hgb hct mcv mch mchc ///
               plt mpv pdw pct

foreach v of local biomvars {
    capture confirm variable `v'
    if _rc==0 {
        destring `v', replace force
        * 负值设缺失
        replace `v' = . if `v' < 0
        * 常见"错误码"转缺失
        replace `v' = . if inlist(`v', 7777, 8888, 9999, 777, 888, 999)
    }
}
* 统一四分类标签：1=高于正常, 2=低于正常, 3=正常
* 说明：缺失值保持 .

* ---- 白蛋白 35–55 g/L
gen alb_cat = .
replace alb_cat = 1 if alb > 55
replace alb_cat = 2 if alb < 35
replace alb_cat = 3 if inrange(alb,35,55)
label values alb_cat triocat
label var alb_cat "白蛋白分类(高/低/正常)"

* ---- 血糖 3.9–6.1 mmol/L
gen glu_cat = .
replace glu_cat = 1 if glu > 6.1
replace glu_cat = 2 if glu < 3.9
replace glu_cat = 3 if inrange(glu,3.9,6.1)
label values glu_cat triocat
label var glu_cat "血糖分类(高/低/正常)"

* ---- 血尿素氮 2.5–6.3 mmol/L
gen bun_cat = .
replace bun_cat = 1 if bun > 6.3
replace bun_cat = 2 if bun < 2.5
replace bun_cat = 3 if inrange(bun,2.5,6.3)
label values bun_cat triocat
label var bun_cat "血尿素氮分类(高/低/正常)"

* ---- 肌酐 50–120 mmol/L
gen crea_cat = .
replace crea_cat = 1 if crea > 120
replace crea_cat = 2 if crea < 50
replace crea_cat = 3 if inrange(crea,50,120)
label values crea_cat triocat
label var crea_cat "肌酐分类(高/低/正常)"

* ---- 总胆固醇 0–5.18 mmol/L
gen cho_cat = .
replace cho_cat = 1 if cho > 5.18
replace cho_cat = 2 if cho < 0
replace cho_cat = 3 if inrange(cho,0,5.18)
label values cho_cat triocat
label var cho_cat "胆固醇分类(高/低/正常)"

* ---- 甘油三酯 0.57–1.71 mmol/L
gen tg_cat = .
replace tg_cat = 1 if tg > 1.71
replace tg_cat = 2 if tg < 0.57
replace tg_cat = 3 if inrange(tg,0.57,1.71)
label values tg_cat triocat
label var tg_cat "甘油三酯分类(高/低/正常)"

* ---- 糖化血清蛋白 1.08–2.10 mmol/L
gen gsp_cat = .
replace gsp_cat = 1 if gsp > 2.10
replace gsp_cat = 2 if gsp < 1.08
replace gsp_cat = 3 if inrange(gsp,1.08,2.10)
label values gsp_cat triocat
label var gsp_cat "糖化血清蛋白分类(高/低/正常)"

* ---- 尿酸 180–440 umol/L
gen ua_cat = .
replace ua_cat = 1 if ua > 440
replace ua_cat = 2 if ua < 180
replace ua_cat = 3 if inrange(ua,180,440)
label values ua_cat triocat
label var ua_cat "尿酸分类(高/低/正常)"

* ---- HDL-C 1.03–1.55 mmol/L（高密度脂蛋白）
gen hdlc_cat = .
replace hdlc_cat = 1 if hdlc > 1.55
replace hdlc_cat = 2 if hdlc < 1.03
replace hdlc_cat = 3 if inrange(hdlc,1.03,1.55)
label values hdlc_cat triocat
label var hdlc_cat "HDL-C分类(高/低/正常)"

* ---- SOD 42.5–76.5 IU/ml
gen sod_cat = .
replace sod_cat = 1 if sod > 76.5
replace sod_cat = 2 if sod < 42.5
replace sod_cat = 3 if inrange(sod,42.5,76.5)
label values sod_cat triocat
label var sod_cat "SOD分类(高/低/正常)"

* ---- WBC 4–10 (10^9/L)
gen wbc_cat = .
replace wbc_cat = 1 if wbc > 10
replace wbc_cat = 2 if wbc < 4
replace wbc_cat = 3 if inrange(wbc,4,10)
label values wbc_cat triocat
label var wbc_cat "白细胞计数分类(高/低/正常)"

* ---- 淋巴细胞计数 0.8–4.0 (10^9/L)
gen lymph_cat = .
replace lymph_cat = 1 if lymph > 4
replace lymph_cat = 2 if lymph < 0.8
replace lymph_cat = 3 if inrange(lymph,0.8,4)
label values lymph_cat triocat
label var lymph_cat "淋巴细胞计数分类(高/低/正常)"

* ---- 淋巴细胞百分比 17–50 (%)
gen lymph_a_cat = .
replace lymph_a_cat = 1 if lymph_a > 50
replace lymph_a_cat = 2 if lymph_a < 17
replace lymph_a_cat = 3 if inrange(lymph_a,17,50)
label values lymph_a_cat triocat
label var lymph_a_cat "淋巴细胞百分比分类(高/低/正常)"

* ---- RBC 3.5–5.5 (10^12/L)
gen rbc_cat = .
replace rbc_cat = 1 if rbc > 5.5
replace rbc_cat = 2 if rbc < 3.5
replace rbc_cat = 3 if inrange(rbc,3.5,5.5)
label values rbc_cat triocat
label var rbc_cat "红细胞计数分类(高/低/正常)"

* ---- HGB 110–160 g/L
gen hgb_cat = .
replace hgb_cat = 1 if hgb > 160
replace hgb_cat = 2 if hgb < 110
replace hgb_cat = 3 if inrange(hgb,110,160)
label values hgb_cat triocat
label var hgb_cat "血红蛋白分类(高/低/正常)"

* ---- HCT 37–50 (%)
gen hct_cat = .
replace hct_cat = 1 if hct > 50
replace hct_cat = 2 if hct < 37
replace hct_cat = 3 if inrange(hct,37,50)
label values hct_cat triocat
label var hct_cat "红细胞压积分(高/低/正常)"

* ---- MCV 80–100 fL
gen mcv_cat = .
replace mcv_cat = 1 if mcv > 100
replace mcv_cat = 2 if mcv < 80
replace mcv_cat = 3 if inrange(mcv,80,100)
label values mcv_cat triocat
label var mcv_cat "平均红细胞体积分(高/低/正常)"

* ---- MCH 26–38 pg
gen mch_cat = .
replace mch_cat = 1 if mch > 38
replace mch_cat = 2 if mch < 26
replace mch_cat = 3 if inrange(mch,26,38)
label values mch_cat triocat
label var mch_cat "平均红细胞血红蛋白量分(高/低/正常)"

* ---- MCHC 300–360 g/L
gen mchc_cat = .
replace mchc_cat = 1 if mchc > 360
replace mchc_cat = 2 if mchc < 300
replace mchc_cat = 3 if inrange(mchc,300,360)
label values mchc_cat triocat
label var mchc_cat "平均红细胞血红蛋白浓度分(高/低/正常)"

* ---- 血小板 PLT 100–300 (10^9/L)
gen plt_cat = .
replace plt_cat = 1 if plt > 300
replace plt_cat = 2 if plt < 100
replace plt_cat = 3 if inrange(plt,100,300)
label values plt_cat triocat
label var plt_cat "血小板计数分类(高/低/正常)"

* ---- MPV 7–13 fL
gen mpv_cat = .
replace mpv_cat = 1 if mpv > 13
replace mpv_cat = 2 if mpv < 7
replace mpv_cat = 3 if inrange(mpv,7,13)
label values mpv_cat triocat
label var mpv_cat "平均血小板体积分(高/低/正常)"

* ---- PDW 10–18 fL
gen pdw_cat = .
replace pdw_cat = 1 if pdw > 18
replace pdw_cat = 2 if pdw < 10
replace pdw_cat = 3 if inrange(pdw,10,18)
label values pdw_cat triocat
label var pdw_cat "血小板体积分布宽度分(高/低/正常)"

* ---- PCT 0.10–0.35 %
gen pct_cat = .
replace pct_cat = 1 if pct > 0.35
replace pct_cat = 2 if pct < 0.10
replace pct_cat = 3 if inrange(pct,0.10,0.35)
label values pct_cat triocat
label var pct_cat "血小板压积分(高/低/正常)"

* 对偏态变量做温和 winsor（P1/P99）
local skewvars glu tg gsp crphs ua mda ualb ucr ualbbyucr plt mpv pdw
foreach v of local skewvars {
    capture noisily su `v', detail
    if _rc==0 & r(N)>0 {
        local p1 = r(p1)
        local p99 = r(p99)
        replace `v' = `p1'  if `v' < `p1'
        replace `v' = `p99' if `v' > `p99'
    }
}

* 标准化 Z 分数（便于回归）
foreach v of local biomvars {
    capture confirm variable `v'
    if _rc==0 {
        egen `v'_z = std(`v')
    }
}

* 派生指标
capture confirm variable cho
capture confirm variable hdlc
if !_rc {
    gen non_hdl = cho - hdlc if !missing(cho, hdlc)
    label var non_hdl "Non-HDL cholesterol"
    gen ln_non_hdl = ln(non_hdl) if non_hdl>0
}
capture confirm variable tg
if !_rc & "`: list hdlc in varlist'" != "" {
    gen aip = ln(tg/hdlc) if tg>0 & hdlc>0
    label var aip "AIP ln(TG/HDL-C)"
}

* 变量标签（按你的中文标签补齐）
label var alb       "白蛋白"
label var glu       "血糖"
label var bun       "血尿素氮"
label var crea      "肌酐"
label var cho       "胆固醇"
label var tg        "甘油三酯"
label var gsp       "糖化血清蛋白"
label var crphs     "超敏C-反应蛋白"
label var ua        "尿酸"
label var hdlc      "高密度脂蛋白胆固醇"
label var sod       "超氧化物歧化酶(SOD)活性"
label var mda       "丙二醛"
label var vd3       "维生素D3"
label var vitb12    "维生素B12"
label var ualb      "尿微量白蛋白"
label var ucr       "尿肌酐"
label var ualbbyucr "尿清蛋白/尿肌酐"
label var wbc       "白细胞计数"
label var lymph     "淋巴细胞计数"
label var lymph_a   "淋巴细胞百分比"
label var rbc       "红细胞计数"
label var hgb       "血红蛋白浓度"
label var hct       "红细胞压积"
label var mcv       "平均红细胞体积"
label var mch       "平均红细胞血红蛋白含量"
label var mchc      "平均红细胞血红蛋白浓度"
label var plt       "血小板计数"
label var mpv       "平均血小板体积"
label var pdw       "血小板体积分布宽度"
label var pct       "血小板压积"
label var non_hdl   "非高密度脂蛋白胆固醇"
label var aip      "动脉粥样硬化指数"

* save
order id, first
save "biomarker_2014_clean.dta", replace
export delimited using "../../Working_data/biomarker_2014_clean.csv", replace

display as text "== 生化清洗完成：已保存 biomarker_2014_clean.dta / .csv =="

