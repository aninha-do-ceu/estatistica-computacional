# carregando o dataset fuel2001.csv

library(dplyr)
library(knitr)
library(ggplot2)
library(patchwork)
library(invgamma)

setwd('C:/projeto/estatistica_computacional')

# ------ Baixando a Base --------- #

df_selic <- read.table('taxa_juros_selic.csv',sep=";",header<-TRUE)
df_ibov <- read.table('ibov.csv',sep=",",header<-TRUE)

# ------ Transformando as colunas --------- #

colnames(df_selic) <- c("data", "selic")
colnames(df_ibov) <- c("indice", "data", "ibov")

df_selic <- df_selic[df_selic$selic != 'BCB-Demab', ]
# Substituindo as vírgulas por pontos e convertendo para 'numeric'
df_selic$selic <- as.numeric(gsub(',', '.', df_selic$selic))
# Convertendo a coluna 'data' para o tipo Date
df_selic$data <- as.Date(df_selic$data, format="%d/%m/%Y")

# Convertendo a coluna 'data' para o tipo Date
df_ibov$data <- as.Date(df_ibov$data, format="%d.%m.%Y")
# Convertendo a coluna 'ibov' para numeric
df_ibov$ibov <- as.numeric(df_ibov$ibov)

# ------ Realizando a junção entre 'ibov' e 'selic' com base na coluna 'data'  --------- #
df_series <- merge(df_ibov, df_selic, by = 'data', all = FALSE)
df_series$ano <- format(df_series$data, "%Y")
df_series$ano <- as.factor(df_series$ano)

# ------ Verificando as primeiras linhas e os tipos  --------- #
head(df_series)

str(df_series)

# ------ Estatísticas das colunas (R) --------- #
summary(df_series)

# ------ Iremos fazer um modelo ajustando uma reta para as variáveis FuelC (resposta) e Drivers (explicativa)  --------- #

# ------ Visualizando por Gráfico de Dispersão  --------- #


plot1 <- ggplot(df_series, aes(x = selic, y = ibov)) +
          geom_point(size = 2, colour=as.factor(df_series$ano)) +  # Aumentando o tamanho dos pontos
          labs(
            title = "Gráfico de Dispersão - IBOVESPA x SELIC",
            y = "Índice IBOVESPA",
            x = "Taxa SELIC"
          ) +
          theme_gray() +
          scale_color_discrete(name = "Ano", labels = unique(df_series$ano)) 

plot1
ggsave("series_tempo.png", plot1, width = 8, height = 6)

# ajustando modelo

x = df_series$ibov
y = df_series$selic

model <- lm(y ~ x)
b0_est <- signif(model$coef[[1]], 5)
b1_est <- signif(model$coef[[2]], 5)

x <- df_series$ibov
y <- df_series$selic

df_series$pred = model$fitted.values

# ic para os coeficientes

sq_total <- sum((df_series$selic - mean(df_series$selic))^2)

# sq_reg tamb?m ? 
s_xy <- function(x,y){
  x_barra <- mean(x)
  y_barra <- mean(y)
  
  return(sum((x - x_barra)*(y - y_barra)))
}
sq_reg <- b1_est*s_xy(x,y)

sq_res <- sq_total - sq_reg

# Estimador de sigma^2: Qm_res

n <- length(x)
qm_res <- sq_res/(n-2)
alfa = 0.05
n <- length(x)
t_quant_critico_1 <- qt(alfa/2, n-2)

dp_b0 <- sqrt(qm_res*((1/n) + (mean(x)^2/s_xy(x,x))))
t_b0 <- b0_est/dp_b0 

dp_b1 <- sqrt(qm_res/s_xy(x,x))
t_b1 <- b1_est/dp_b1

qui_quant_critico_1 <- qchisq(alfa/2, n-2)
qui_quant_critico_2 <- qchisq(1 - (alfa/2), n-2)

sig_min <- sq_res/qui_quant_critico_2
sig_max <- sq_res/qui_quant_critico_1


log_post_eta <- function(eta, y, x, b0, b1, alpha, beta) {
  n <- length(y)
  sigma2 <- exp(2 * eta)
  rss <- sum((y - b0 - b1 * x)^2)
  
  logpost <- -(n + 2 * alpha) * eta -
    (rss / 2 + beta) * exp(-2 * eta)
  
  return(logpost)
}

# Função que realiza o processo de simulação para estimar b0 e b1
gibbs_regression <- function(
    n_iter, x, y, b0, b1, s2, alfa, beta, s2_b, tau
) {
  
  n <- length(x)
  
  b0_atual <- b0
  b1_atual <- b1
  eta_atual <- log(sqrt(s2))  # log(sigma)
  
  a <- alfa
  b <- beta
  
  b0_chain <- numeric(n_iter)
  b1_chain <- numeric(n_iter)
  s2_chain <- numeric(n_iter)
  
  aceitos <- 0
  
  for (i in 1:n_iter) {
    
    s2_atual <- exp(2 * eta_atual)
    
    ## --- Gibbs para b0 ---
    s2_0 <- 1 / (n / s2_atual + 1 / s2_b)
    mu_0 <- s2_0 * sum((y - b1_atual * x) / s2_atual)
    
    b0_atual <- rnorm(1, mu_0, sqrt(s2_0))
    
    ## --- Gibbs para b1 ---
    s2_1 <- 1 / (sum(x^2) / s2_atual + 1 / s2_b)
    mu_1 <- s2_1 * sum(x * (y - b0_atual) / s2_atual)
    
    b1_atual <- rnorm(1, mu_1, sqrt(s2_1))
    
    ## --- MH para log(sigma) ---
    eta_prop <- rnorm(1, eta_atual, sd = tau)
    
    log_acc <- log_post_eta(
      eta_prop, y, x, b0_atual, b1_atual, a, b
    ) -
      log_post_eta(
        eta_atual, y, x, b0_atual, b1_atual, a, b
      )
    
    if (log(runif(1)) < log_acc) {
      eta_atual <- eta_prop
      aceitos <- aceitos + 1
    }
    
    ## --- salvar ---
    b0_chain[i] <- b0_atual
    b1_chain[i] <- b1_atual
    s2_chain[i] <- exp(2 * eta_atual)
  }
  
  list(
    aceitos = aceitos,
    taxa_aceitacao = aceitos / n_iter,
    amostras = data.frame(
      b0 = b0_chain,
      b1 = b1_chain,
      s2 = s2_chain
    )
  )
}

# Chama a função
# priori nao informativa
list_resultado <- gibbs_regression(n_iter = 10000, x = x, y = y, b0=0, b1=0, s2=1000, alfa=0.1, beta=0.1, s2_b=1,tau=0.2)

print('Taxa de aceitação:\n')
print(list_resultado$taxa_aceitacao)

resultado <- list_resultado$amostras
resultado$indice <- seq(1:10000)

resultado <- resultado %>%
  filter(indice > 1000) %>%           
  slice(seq(1, n(), by = 5)) 


plot1 <- ggplot(resultado, aes(x = indice, y = b0)) +
  geom_line(color = "#2C7BB6", linewidth = 0.4) +
  labs(
    title = expression("Traço da Cadeia - Coeficiente " * beta[0]),
    x = "Iteração",
    y = expression(beta[0])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

plot2 <- ggplot(resultado, aes(x = indice, y = b1)) +
  geom_line(color = "#2C7BB6", linewidth = 0.4) +
  labs(
    title = expression("Traço da Cadeia - Coeficiente " * beta[1]),
    x = "Iteração",
    y = expression(beta[1])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

plot3 <- ggplot(resultado, aes(x = indice, y = s2)) +
  geom_line(color = "#2C7BB6", linewidth = 0.4) +
  labs(
    title = expression("Traço da Cadeia - Coeficiente " * sigma),
    x = "Iteração",
    y = expression(sigma)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

(plot1 + plot2 + plot3)

ggsave("traco_cadeias.png", (plot1 + plot2 + plot3), width = 12, height = 6)


library(reshape2)

# calcula ACF manualmente
acf_b0 <- acf(resultado$b0, plot=FALSE, lag.max=50)
acf_df <- data.frame(lag = seq(1:51), acf = acf_b0$acf)

plot1 <- ggplot(acf_df, aes(x=lag, y=acf)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_hline(yintercept=c(0, 2/sqrt(nrow(df)), -2/sqrt(nrow(df))),
             linetype="dashed", color="red") +
  labs(title=expression("ACF - Coeficiente " * beta[0]), x="Lag", y="Autocorrelação") +
  theme_minimal()

# calcula ACF manualmente
acf_b1 <- acf(resultado$b1, plot=FALSE, lag.max=50)
acf_df <- data.frame(lag = seq(1:51), acf = acf_b1$acf)

plot2 <- ggplot(acf_df, aes(x=lag, y=acf)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_hline(yintercept=c(0, 2/sqrt(nrow(df)), -2/sqrt(nrow(df))),
             linetype="dashed", color="red") +
  labs(title=expression("ACF - Coeficiente " * beta[1]), x="Lag", y="Autocorrelação") +
  theme_minimal()

# calcula ACF manualmente
acf_s2 <- acf(resultado$s2, plot=FALSE, lag.max=50)
acf_df <- data.frame(lag = seq(1:51), acf = acf_s2$acf)
plot3 <- ggplot(acf_df, aes(x=lag, y=acf)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_hline(yintercept=c(0, 2/sqrt(nrow(df)), -2/sqrt(nrow(df))),
             linetype="dashed", color="red") +
  labs(title=expression("ACF - Coeficiente " * beta[0]), x="Lag", y="Autocorrelação") +
  theme_minimal()

(plot1 + plot2 + plot3)
ggsave("autocorr.png", (plot1 + plot2 + plot3), width = 12, height = 6)

# priori informativa
list_resultado_inf <- gibbs_regression(n_iter = 10000, x = x, y = y, b0=0, b1=0, s2=2, alfa=0.1, beta=0.1, s2_b=1,tau=0.1)

print('Taxa de aceitação:\n')
print(list_resultado$taxa_aceitacao)

resultado_inf <- list_resultado_inf$amostras
resultado_inf$indice <- seq(1:10000)

resultado_inf <- resultado_inf %>%
  filter(indice > 1000) %>%           
  slice(seq(1, n(), by = 5)) 


plot1 <- ggplot(resultado_inf, aes(x = indice, y = b0)) +
  geom_line(color = "#2C7BB6", linewidth = 0.4) +
  labs(
    title = expression("Traço da Cadeia - Coeficiente " * beta[0]),
    x = "Iteração",
    y = expression(beta[0])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

plot2 <- ggplot(resultado_inf, aes(x = indice, y = b1)) +
  geom_line(color = "#2C7BB6", linewidth = 0.4) +
  labs(
    title = expression("Traço da Cadeia - Coeficiente " * beta[1]),
    x = "Iteração",
    y = expression(beta[1])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

plot3 <- ggplot(resultado_inf, aes(x = indice, y = s2)) +
  geom_line(color = "#2C7BB6", linewidth = 0.4) +
  labs(
    title = expression("Traço da Cadeia - Coeficiente " * sigma),
    x = "Iteração",
    y = expression(sigma)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

(plot1 + plot2 + plot3)

ggsave("traco_cadeias_inf.png", (plot1 + plot2 + plot3), width = 12, height = 6)


library(reshape2)

# calcula ACF manualmente
acf_b0 <- acf(resultado_inf$b0, plot=FALSE, lag.max=50)
acf_df <- data.frame(lag = seq(1:51), acf = acf_b0$acf)

plot1 <- ggplot(acf_df, aes(x=lag, y=acf)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_hline(yintercept=c(0, 2/sqrt(nrow(df)), -2/sqrt(nrow(df))),
             linetype="dashed", color="red") +
  labs(title=expression("ACF - Coeficiente " * beta[0]), x="Lag", y="Autocorrelação") +
  theme_minimal()

# calcula ACF manualmente
acf_b1 <- acf(resultado_inf$b1, plot=FALSE, lag.max=50)
acf_df <- data.frame(lag = seq(1:51), acf = acf_b1$acf)

plot2 <- ggplot(acf_df, aes(x=lag, y=acf)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_hline(yintercept=c(0, 2/sqrt(nrow(df)), -2/sqrt(nrow(df))),
             linetype="dashed", color="red") +
  labs(title=expression("ACF - Coeficiente " * beta[1]), x="Lag", y="Autocorrelação") +
  theme_minimal()

# calcula ACF manualmente
acf_s2 <- acf(resultado_inf$s2, plot=FALSE, lag.max=50)
acf_df <- data.frame(lag = seq(1:51), acf = acf_s2$acf)
plot3 <- ggplot(acf_df, aes(x=lag, y=acf)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_hline(yintercept=c(0, 2/sqrt(nrow(df)), -2/sqrt(nrow(df))),
             linetype="dashed", color="red") +
  labs(title=expression("ACF - Coeficiente " * sigma), x="Lag", y="Autocorrelação") +
  theme_minimal()

(plot1 + plot2 + plot3)
ggsave("autocorr_inf.png", (plot1 + plot2 + plot3), width = 12, height = 6)


df_b0 <- rbind(data.frame(b0 = list_resultado$amostras$b0, priori =rep("Fracamente Informativa",10000)),
data.frame(b0 =list_resultado_inf$amostras$b0, priori = rep("Informativa",10000)))

df_b1 <- rbind(data.frame(b1 = list_resultado$amostras$b1, priori =rep("Fracamente Informativa",10000)),
               data.frame(b1 =list_resultado_inf$amostras$b1, priori = rep("Informativa",10000)))
df_s2 <- rbind(data.frame(s2 = list_resultado$amostras$s2, priori =rep("Fracamente Informativa",10000)),
               data.frame(s2 =list_resultado_inf$amostras$s2, priori = rep("Informativa",10000)))


plot1 <-ggplot(data =resultado, aes(x=b0)) +
  geom_density(alpha=0.4) +
  geom_vline(xintercept=b0_est, color="red", linetype="dashed", size=1) +
  labs(title=expression("Densidade Marginal de " * beta[0] * " - Priori Informativa Fraca"),
       x="Valor",
       y="Densidade") +
  xlim(0,0.2)+
  theme_minimal() 

plot2 <-ggplot(data =resultado_inf, aes(x=b0)) +
  geom_density(alpha=0.4) +
  geom_vline(xintercept=b0_est, color="red", linetype="dashed", size=1) +
  labs(title=expression("Densidade Marginal de " * beta[0] * " - Priori Informativa"),
       x="Valor",
       y="Densidade") +
  xlim(0,0.2)+
  theme_minimal() 

plot1/plot2

ggsave("dens_b0.png", plot1/plot2)

plot1 <-ggplot(data =resultado, aes(x=b1)) +
  geom_density(alpha=0.4) +
  geom_vline(xintercept=b1_est, color="red", linetype="dashed", size=1) +
  labs(title=expression("Densidade Marginal de " * beta[1] * " - Priori Informativa Fraca"),
       x="Valor",
       y="Densidade") +
  xlim(-0.0005,0.0005)+
  theme_minimal() 

plot2 <-ggplot(data =resultado_inf, aes(x=b1)) +
  geom_density(alpha=0.4) +
  geom_vline(xintercept=b1_est, color="red", linetype="dashed", size=1) +
  labs(title=expression("Densidade Marginal de " * beta[1] * " - Priori Informativa"),
       x="Valor",
       y="Densidade") +
  xlim(-0.0005,0.0005)+
  theme_minimal() 

plot1/plot2

ggsave("dens_b1.png", plot1/plot2)

plot1 <-ggplot(data =resultado, aes(x=s2)) +
  geom_density(alpha=0.4) +
  labs(title=expression("Densidade Marginal de " * sigma * " - Priori Informativa Fraca"),
       x="Valor",
       y="Densidade") +
  xlim(0,0.002) +
  theme_minimal() 

plot2 <-ggplot(data =resultado_inf, aes(x=s2)) +
  geom_density(alpha=0.4) +
  labs(title=expression("Densidade Marginal de " * sigma * " - Priori Informativa"),
       x="Valor",
       y="Densidade") +
  xlim(0,0.002) +
  theme_minimal() 

plot1/plot2

ggsave("dens_s2.png", plot1/plot2)



log_lik_pointwise <- function(y, x, b0_chain, b1_chain, s2_chain) {
  n_iter <- length(b0_chain)
  n <- length(y)
  
  # matriz de log-likelihood: linhas = iterações, colunas = observações
  log_lik_mat <- matrix(NA, nrow=n_iter, ncol=n)
  
  for (s in 1:n_iter) {
    mu <- b0_chain[s] + b1_chain[s] * x
    sigma <- sqrt(s2_chain[s])
    log_lik_mat[s, ] <- dnorm(y, mean=mu, sd=sigma, log=TRUE)
  }
  
  return(log_lik_mat)
}


compute_WAIC <- function(log_lik_mat) {
  # lppd
  lppd <- sum(log(colMeans(exp(log_lik_mat))))
  
  # p_WAIC
  p_WAIC <- sum(apply(log_lik_mat, 2, var))
  
  # WAIC
  WAIC <- -2 * (lppd - p_WAIC)
  
  return(list(WAIC=WAIC, lppd=lppd, p_WAIC=p_WAIC))
}

# matriz de log-likelihood
log_lik <- log_lik_pointwise(
  y=y, x=x,
  b0_chain=resultado$b0,
  b1_chain=resultado$b1,
  s2_chain=resultado$s2
)

# WAIC
waic <- compute_WAIC(log_lik)
print(waic$WAIC)



log_lik <- log_lik_pointwise(
  y=y, x=x,
  b0_chain=resultado_inf$b0,
  b1_chain=resultado_inf$b1,
  s2_chain=resultado_inf$s2
)


# WAIC
waic <- compute_WAIC(log_lik)
print(waic$WAIC)


# Função para calcular a moda aproximada usando KDE
calc_moda <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

# Função para calcular intervalo de 95%
calc_ic95 <- function(x) {
  quantile(x, probs=c(0.025, 0.975))
}
summarize_chain <- function(chain, name, priori, waic=0) {
  ic <- calc_ic95(chain)
  
  data.frame(
    priori = priori,
    parametro = name,
    media = mean(chain),
    mediana = median(chain),
    moda = calc_moda(chain),
    ic_95_inf = ic[1],
    ic_95_sup = ic[2],
    WAIC = waic
  )
}
# Cadeia 1
b0_chain1 <- resultado_inf$b0
b1_chain1 <- resultado_inf$b1
s2_chain1 <- resultado_inf$s2

# Cadeia 2
b0_chain2 <- resultado$b0
b1_chain2 <- resultado$b1
s2_chain2 <- resultado$s2

# Resumo para primeira cadeia
df1 <- rbind(
  summarize_chain(b0_chain1, "b0", priori= "Informativa", waic=-817.8362),
  summarize_chain(b1_chain1, "b1", priori="Informativa", waic=-817.8362),
  summarize_chain(s2_chain1, "s2", priori="Informativa", waic=-817.8362)
)

# Resumo para segunda cadeia
df2 <- rbind(
  summarize_chain(b0_chain2, "b0", priori= "Não Informativa", waic=-817.8869),
  summarize_chain(b1_chain2, "b1", priori= "Não Informativa", waic=-817.8869),
  summarize_chain(s2_chain2, "s2", priori= "Não Informativa", waic=-817.8869)
)

# Adicionando identificação da cadeia
df1$cadeia <- "resultado_inf"
df2$cadeia <- "resultado"

# Juntar os dois data frames
df_final <- rbind(df1, df2)

rbind(
  df_final,
  df3
)

df3 <- data.frame(priori = c("Modelo Original","Modelo Original","Modelo Original"),
                  parametro = c("b0","b1","s2"),
                  media = c(b0_est, b1_est, "-"),
                  mediana = c("-","-","-"),
                  moda = c("-","-","-"),
                  ic_95_inf = c(t_b0 - dp_b0,t_b1 - dp_b1, sig_min),
                  ic_95_sup = c(t_b0 + dp_b0,t_b1 + dp_b1, sig_max),
                  WAIC=c("-","-","-"),
                  cadeia = c("res","",""))
