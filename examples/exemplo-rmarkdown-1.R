# ----------------------------------
# Porto Alegre, 15 de outubro de 2018
# Script de exemplo Rmarkdown
# Autor: Rodrigo

# ----------------------------------
# Carregando pacotes

library(gapminder)
library(dplyr)
library(readr)
library(ggplot2)
library(compareGroups)

# ----------------------------------
# Manipulação de dados

gapminder <- gapminder %>%
  mutate(pop_m = pop/1e6)

gapminder$gdpPercap.cat <- cut(gapminder$gdpPercap,
                               breaks = c(0, 1005, 3955, 12235, Inf),
                               labels = c("Baixa-renda", "Renda média-baixa",
                                          "Renda média-alta", "Renda alta"))

gapminder07 <- gapminder %>%
  filter(year == 2007)

# ----------------------------------
# Uma tabela descritiva

# summary(gapminder07)

gapm.df <- as.data.frame(gapminder07)

res <- compareGroups(continent ~ . - year - country, data = gapm.df)
restab <- createTable(res)
print(restab, which.table = 'descr')

# ----------------------------------
# Um resumo numérico

max(gapminder07$lifeExp)
gapminder07$country[which.max(gapminder07$lifeExp)]

# ----------------------------------
# Um gráfico descritivo

p <- ggplot(data = gapminder07,
            mapping = aes(x = gdpPercap, y = lifeExp,
                          color = continent, size = pop_m)) + 
  geom_point() +
  labs(x = "Renda per capita (US$)",
       y = "Expectativa de vida (anos)",
       color = "Continente", size = "População/1 milhão") +
  theme_bw()
p

# ----------------------------------
# Um modelo de regressão linear

mod1 <- lm(lifeExp ~ gdpPercap, data = gapminder07)
summary(mod1)$coef

mod2 <- lm(lifeExp ~ gdpPercap + continent, data = gapminder07)
summary(mod2)$coef

# ----------------------------------
# Um gráfico exploratório

p + geom_smooth(mapping = aes(x = gdpPercap, y = lifeExp, color = NULL),
                method = "loess")

