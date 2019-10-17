# ----------------------------------
# Carregando pacotes

library(knitr)
library(gapminder)
library(dplyr)

as.character(unique(gapminder$year))
  
filelist <- unique(gapminder$year)

gapminder <- gapminder %>%
  mutate(pop_m = pop/1e6)

gapminder$gdpPercap.cat <- cut(gapminder$gdpPercap,
                               breaks = c(0, 1005, 3955, 12235, Inf),
                               labels = c("Baixa-renda", "Renda média-baixa",
                                          "Renda média-alta", "Renda alta"))

# loop through the file list to read in data and clean it up

for (file in filelist) {
  
  fp <- paste("C:\\Users\\Rodrigo\\Documents\\PintandoEBordando\\ArquivosR\\relatorio_gapminder\\", file, sep="")
  # print(fp)

  # read in continents
  gpm <- gapminder %>% 
    filter(year == file)
  
  rmarkdown::render(input = "relatorio_gapminder\\rel_gap_template.Rmd", 
                    output_format = "html_document",
                    output_file = paste0("Relatorio_gapminder_", file, ".html"),
                    output_dir = "relatorio_gapminder\\mult_relatorios",
                    encoding = "UTF-8")
  
}