# ------------- config -------------
START, END = "2004-01", "2023-12"
COUNTRIES = {"CZ": "Czechia", "SK": "Slovakia"}

# eurostat datasets & filters you want
EUROSTAT_SERIES = {
    # name        dataset               filter dict  (see Eurostat API cheatsheet)
    "hicp_idx":   ("prc_hicp_midx",     dict(coicop="CP00", unit="I15")),         # all-items index
    "unemp_sa":   ("une_rt_m",          dict(sex="T", age="Y15-74", seasonaladj="SA")),
    "gov_debt":   ("gov_10q_ggdebt",    dict(unit="PC_GDP", sector="S13")),
    "fdi_in":     ("bop_fdi6_flow",     dict(direction="IN", partner="WLD", assetliab="LIAB")),
    "fdi_out":    ("bop_fdi6_flow",     dict(direction="OUT", partner="WLD", assetliab="LIAB")),
}

ECB_SDW = {
    "policy_rate_SK": "FM.M.U2.EUR.4F.IR.M",      # main refi (euro area proxy for SK)
}

CNB_SERIES = {
    "policy_rate_CZ": "https://www.cnb.cz/en/statistics/monetary_statistics/"  # 2-week repo JSON; endpoint below
                  "financial_markets/fin_data/repo_2w.json"
}