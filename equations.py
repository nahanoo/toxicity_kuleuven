from sympy import symbols, solve, Eq, init_printing
from model import Salmonella, E_coli, Experiment

e_coli = E_coli()
salm = Salmonella()
exp = Experiment()

init_printing()
(
    e,
    s,
    R,
    r,
    u,
    Km,
    KT,
    Y,
    cefo,
    chloram,
    dE,
    dS,
    dR,
    dCefo,
    dChloram,
    ue,
    us,
    f,
    a,
    j,
    e_sus,
    s_sus,
) = symbols(
    "e,s,R,r,u,Km,KT,Y,cefo,chloram,dE,dS,dR,dCefo,dChloram,ue,us,f,a,j,e_sus,s_sus"
)

e_sus = Eq(e_sus, e * r * R / (R + Km) - e * u * cefo / (cefo + KT))
e_min_u = solve(e_sus.rhs, u)[0].subs(
    {"R": exp.M, "r": e_coli.r, "cefo": exp.M_cefo, "Km": e_coli.Km, "KT": e_coli.KT}
)
s_sus = Eq(s_sus, s * r * R / (R + Km) - s * u * chloram / (chloram + KT))
s_min_u = solve(s_sus.rhs, u)[0].subs(
    {"R": exp.M, "r": salm.r, "chloram": exp.M_chloram, "Km": salm.Km, "KT": salm.KT}
)
