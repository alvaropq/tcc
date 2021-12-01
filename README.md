# PROPOSTA DE UM PIPELINE BASEADO EM APRENDIZADO DE MÁQUINA PARA GERAÇÃO E ANÁLISE DE SEQUÊNCIAS BIOLÓGICAS
## Álvaro Pedroso Queiroz - UTFPR
## Danilo Sipoli Sanches - UTFPR
O avanço da pesquisa na área de bioinformática tem sido emergente. O tamanho dos dados acumulados em vários projetos de sequenciamento está aumentando exponencialmente, e assim, técnicas computacionais envolvendo algoritmos de classificação e agrupamento foram propostos para reduzir as dificuldades encontradas em métodos experimentais. Nesse contexto, várias metodologias foram desenvolvidas com base em genômica comparativa para análise e anotações de suas funcionalidades. Dessa forma, este trabalho desenvolveu um _pipeline_ baseado em aprendizado de máquina para a realização de análises de sequências biológicas. O _pipeline_ proposto utiliza 7 técnicas de extração de características diferentes, com descritores matemáticos e convencionais, aplicados em dados de genomas de cepas bacterianas do gênero _bacilus_. Os experimentos realizados apresentaram resultados promissores, com informações de similaridade próximas aos encontrados utilizando outras ferramentas encontradas em trabalhos correlatos. Também é apresentado a importância das características utilizadas nos modelos criados para buscar informações que possam ser de grande importância para biólogos e pesquisadores sobre as diferenças entre os genomas analisados.

### Estrutura do repositório
#### Codes
Apresenta os códigos realizados em notebook com a linguagem python para combinação das bases de dados extraídas pelo script de feature extraction e o codigo do pipeline gerado para a realização do cálculo de distâncias euclidianas, agrupamento hierarquico, classificação e feature importance.
#### Feature Extraction
Apresenta um script em linguagem python que utiliza o código fonte de técnicas de extração de caracterśticas obtidas através da ferramenta MathFeature
#### Results
Apresenta os resultados obtidos utilizando as ferramentas SplitsTree4, Gegenees e o pipeline criado, utilizando os genomas da pasta Sequences.
#### Sequences
Apresenta os genomas das cepas bacterianas utilizadas nos experimentos.
