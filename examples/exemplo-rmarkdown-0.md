# Um documento em Markdown

## Sobre o Markdown

O Markdown é uma linguagem de marcação muito simples, desenvolvida por
John Gruber.

A ideia básica por trás da linguagem é fazer com que o escritor se
preocupe mais com o **conteúdo** do texto do que com a *formatação*.

## Mais um título

Aqui vamos tentar descrever uma análise.

## Simulando variáveis aleatórias

No R podemos simular valores de uma distribuição normal padrão através
da função `rnorm()`.

Seja $X \sim \text{N}(0,1)$, então para gerar 30 valores dessa variável
aleatório normal, fazemos

```
(x <- rnorm(30))
```

Com o resultado dessa simulação, podemos calcular a média e a variância
dessa VA $X$ para conferir se o resultado fica próximo de 0 e 1,
respectivamente.

Também podemos fazer um histograma dessa VA $X$ simulada

```
hist(x)
```