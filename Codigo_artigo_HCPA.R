# Código-fonte para a revista do HCPA
# Estimativas de UTI/COVID-19 atualizadas até 03/08/20
#https://www.ufrgs.br/covidpoa/?p=672

#### CSV PMPA UTI
# https://infografico-covid.procempa.com.br

require(deSolve)
require(ggplot2)
library(lubridate)
require(grid)
options(stringsAsFactors = FALSE)
###########
arq_UTI_pmpa = "1_UTI.csv"
arq_UTI_pmpa_ate_31_12 = "1_UTI_31_12.csv"
arq_Hosp_pmpa = "2_Hosp.csv"

caminho = "Dados_PMPA/"

UTI_pmpa = read.csv(file = paste(caminho, arq_UTI_pmpa, sep = ""))
UTI_pmpa_ate_31_12 = read.csv(file = paste(caminho, arq_UTI_pmpa_ate_31_12, sep = ""))
hosp_pmpa = read.csv(file = paste(caminho, arq_Hosp_pmpa, sep = ""))

index_iqual_days_pmpa = which(UTI_pmpa[,1] %in% hosp_pmpa[,1])

number_reported_full = UTI_pmpa[index_iqual_days_pmpa,2] + hosp_pmpa[,2]
number_reported = number_reported_full[-c(1:24)]
time = c(1:(length(number_reported)))
reported_data = data.frame(time, number_reported)

### Leitos
total_leitos_covid = 383
prop_UTI_H = 0.415

#### Datas
date_ini = dmy("25/Abr/2020") 
num_days = 250 # mais o primeiro dia
series_days = date_ini + days(0:num_days)
times  = seq(0, num_days, by = 1)

num_days_hospi = length(number_reported)
date_ini_hospi = dmy("25/Abr/2020") 

# Parâmetros
sigma  = 1/5.1
gammah = 1/4.86 
gammad = 0 
gammar = 1/5.6
etad   = 1/12.98 
etar   = 1/10.22  
Theta  = 0.01 
Lambda = 0.175  

N = 1014009
I = N * 0.0017 
E = 1.5 * I 
H = 49
D = 12
R = 0
S = N - I - H - E - D - R

initial_state_values = c(S = S, E = E, I = I, H = H, D = D, R = R, incid_acum_I = I, exit_I = 0)

times = seq(from = 0, to = num_days, by = 1)

# SEIHDR 
seihdr_model = function(times, state, parameters) {  
  
  with(as.list(c(state, parameters)), {
    
    incid_I = (sigma * E)
    exit_I  = (Theta * gammah * I) + 
              ((1 - Theta) * (1 - Lambda) * gammar * I) + 
              ((1 - Theta) * Lambda *  gammad * I)
    
    dS = - (beta * I) * (S / N)
    
    dE =   (beta * I) * (S / N) - (sigma * E)
    
    dI =   incid_I - exit_I
    
    dH =   (Theta * gammah * I) - (Lambda * etad * H) - 
      ((1 - Lambda) * etar * H)
    
    dD =   ((1 - Theta) * Lambda *  gammad * I) + 
      (Lambda * etad * H) 
    
    dR =   ((1 - Theta) * (1 - Lambda) * gammar * I) + 
      ((1 - Lambda) * etar * H)
    
    return(list(c(dS, dE, dI, dH, dD, dR, incid_I, exit_I))) 
  })
}

# Função Distância
loglikelihood_fun = function(parameters, dat) { 
  
  beta = parameters[1]
  
  output = as.data.frame(ode(y = initial_state_values, 
                             times = times, 
                             func = seihdr_model,
                             parms = c(beta = beta)))  
  
  # Calcula log-likelihood 
  LL = -sum(dpois(x = dat$number_reported, lambda = output$H[output$time %in% dat$time], log = TRUE))
  
  return(LL) 
}

# Otimização
param_otimim = optim(par = c(0.4),
                     fn = loglikelihood_fun,
                     dat = reported_data, 
                     control = list(fnscale = 1),
                     method = 'Brent',
                     lower = 0,
                     upper = 2) 

beta_estimated = round(param_otimim$par[1], 4)
######################
parameters = c(beta   = beta_estimated, 
               sigma  = sigma,
               gammah = gammah, 
               gammad = gammad, 
               gammar = gammar,
               etad   = etad,
               etar   = etar, 
               Theta  = Theta,
               Lambda = Lambda)

output = as.data.frame(ode( y     = initial_state_values, 
                            times = times, 
                            func  = seihdr_model,
                            parms = parameters))

# Nova coluna de Incidentes I
output$incid_I = c(0,diff(output$incid_acum_I, lag = 1))

### Datas
# Troca o indice do tempo
output$time = series_days
time = date_ini_hospi + days(1:num_days_hospi - 1)
reported_data = data.frame(time, number_reported)

last_day_hosp = tail(time, n = 1)
d = day(last_day_hosp)
m = months(last_day_hosp, abbreviate = T)
y = year(last_day_hosp)
day_month_year_last_day_hosp = paste(d,"/",m,"/",y)

UTI_df = data.frame(UTI_pmpa[-c(1:37),])   
date_ini_UTI_PMPA = dmy("25/Abr/2020")
num_days_UTI_pmpa = length(UTI_df$category)
UTI_days_PMPA = date_ini_UTI_PMPA + days(1:num_days_UTI_pmpa - 1)
UTI_df$category = UTI_days_PMPA
UTI_SEIHDR_prop_UTI_H = output$H * prop_UTI_H
UTI_SEIHDR = data.frame(series_days,UTI_SEIHDR_prop_UTI_H)

### UTI até 31/12
UTI_df_ate_31_12 = data.frame(UTI_pmpa_ate_31_12[-c(1:37),])   
num_days_UTI_pmpa_ate_31_12 = length(UTI_df_ate_31_12$category)
UTI_days_PMPA_ate_31_12 = date_ini_UTI_PMPA + days(1:num_days_UTI_pmpa_ate_31_12 - 1)
UTI_df_ate_31_12$category = UTI_days_PMPA_ate_31_12
UTI_pmpa_ate_31_12_data = data.frame(series_days, UTI_df_ate_31_12$Leitos.UTI)

###
base_UTI = ggplot() +
  geom_line(data = UTI_SEIHDR, aes(x = series_days, y = UTI_SEIHDR_prop_UTI_H), colour = '#56B4E9', size = 2) +
  geom_point(data = UTI_pmpa_ate_31_12_data, aes(x = series_days, y = UTI_df_ate_31_12.Leitos.UTI, colour = "Depois do dia 24/08"), size = 2) +
  scale_color_manual(values = c("black", "red"))+
  geom_point(data = UTI_df, aes(x = category, y = Leitos.UTI, colour = "Até o dia 24/08"), size = 2) +
  xlab("Número observado de pacientes em UTIs/COVID-19") +
  ylab("Número previsto de INTERNADOS por COVID-19 em UTIs") +                                 
  labs(title = paste("Modelo Calibrado até",day_month_year_last_day_hosp), colour = "") +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(angle = 30, vjust = 0.5, size = 13)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(angle = 0, vjust = 0.5, size = 13)) +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 18))+
  geom_hline(yintercept = total_leitos_covid, colour = "orange", size = 1.5, linetype = "dashed") +
  annotate("text", x = dmy("10/jun/2020"), y = total_leitos_covid + 15, label =  paste("Máximo de leitos UTI Porto Alegre/RS,",total_leitos_covid), colour = "orange")+
  geom_text()

# Salva em SVG 
#svg("fig-uti-3-08.svg",width = 10, height = 5)
base_UTI + scale_x_date(date_breaks = '2 week', date_labels = "%d %b")
#dev.off() 
#############
# R0 (next-generation matrix)

F1 = quote(beta * S * I / N)
F2 = 0
F3 = 0

###################################################
Vm1 = quote(sigma * E)
Vm2 = quote(Theta * gammah * I + (1 - Theta) * 
              (1 - Lambda) * gammar * I + (1 - Theta) * 
              Lambda * gammad * I)
Vm3 = quote(Lambda * etad * H + (1 - Lambda) * etar * H)
###################################################
Vp1 = 0
Vp2 = quote(sigma * E)
Vp3 = quote(Theta * gammah * I)

###################################################
V1 = substitute(a - b, list(a = Vm1, b = Vp1))
V2 = substitute(a - b, list(a = Vm2, b = Vp2))
V3 = substitute(a - b, list(a = Vm3, b = Vp3))
###################################################
f11 = D(F1, "E"); f12 = D(F1, "I"); f13 = D(F1, "H")
f21 = D(F2, "E"); f22 = D(F2, "I"); f23 = D(F2, "H")
f31 = D(F3, "E"); f32 = D(F3, "I"); f33 = D(F3, "H")
v11 = D(V1, "E"); v12 = D(V1, "I"); v13 = D(V1, "H")
v21 = D(V2, "E"); v22 = D(V2, "I"); v23 = D(V2, "H")
v31 = D(V3, "E"); v32 = D(V3, "I"); v33 = D(V3, "H")
###################################################
paras = list(N = N, S = N, E = 0, I = 0, 
             H = 0, D = 0, R = 0, 
             beta = beta_estimated,
             sigma = sigma, Theta = Theta, 
             Lambda = Lambda, gammah = gammah,
             gammad = gammad, gammar = gammar, 
             etad = etad, etar = etar)

f = with(paras, 
         matrix(c(eval(f11), eval(f12), eval(f13),
                  eval(f21), eval(f22), eval(f23),
                  eval(f31), eval(f32), eval(f33)),
                nrow = 3, byrow = T))

v = with(paras, 
         matrix(c(eval(v11), eval(v12), eval(v13),
                  eval(v21), eval(v22), eval(v23),
                  eval(v31), eval(v32), eval(v33)),
                nrow = 3, byrow = T))
###################################################
R0 = max(eigen(f %*% solve(v))$values)
R0
###################################################
# R effective
Re = R0 * (output$S[output$time == last_day_hosp]/N)
Re
beta_estimated
#########
i = num_days + 1

paste('pop total Ini = ', N)
paste('pop total Fim = ', round(output$S[i] + output$E[i] + output$I[i] + output$H[i] + output$D[i] + output$R[i]))
paste('H max = ', round(max(output$H), 2))
paste('dia = ', output[which.max(output$H),1])
paste('UTI max = ', round(max(UTI_SEIHDR$UTI_SEIHDR_prop_UTI_H), 2))
paste('I max = ', round(max(output$I), 2))
paste('dia = ', output[which.max(output$I),1])
