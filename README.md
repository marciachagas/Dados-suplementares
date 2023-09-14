# Dados-suplementares

Neste repositório estão os arquivos da dissertação "**Análises de Fenótipos Através de Bacias de Atração do Modelo Booleano da Rede de Regulação Gênica da _Pseudomonas aeruginosa_ CCBH4851**" (Chagas, 2023). São eles:

-	O arquivo da rede completa da CCBH-2022 em CSV;
-	Todos os códigos-fonte:
    -	o	da análise estrutural em R e para encontrar RBH em python, implementados por Medeiros *et al*.(13) (GNR.R e RBH.py, respectivamente)), 
    -	da binarização dos dados de RNA-seq em R (bin.R),
    -	para construção do modelo booleano em python (GRN_to_ASSA_0.2.py),
    -	da identificação das bacias de atração em python (find_basin.py);
-	Arquivo com instruções detalhadas e linhas de comando para construção do modelo booleano e para a simulação da trajetória, por Marcelo Trindade dos Santos (commands.txt);
-	Dados de RNA-seq normalizados, por Felicita Mabel (norm.zip);
-	Dados binarizados dos genes da sub-rede núcleo (genes-bin.zip);
-	Os resultados da simulação de trajetória e as bacias de atração (sim.zip);
-	Tabelas adicionais e figuras principais com maior resolução (figs.zip). Depois do download, dar zoom nas imagens para ver os detalhes.

## **Informações sobre os experimentos de RNA-seq _bulk_:**
Para as curvas de crescimento com antibióticos:
-    A P.aeruginosa foi incubada por 10h (curva de crescimento) em meio Mueller Hinton. Após esse período, foi adicionado os antimicrobianos:
    -    Imipenem = 128 mg/l,
    -    Polimixina B = 1 mg/ml,
    -    Controle (sem adição de antibiótico);

-    As culturas foram incubadas por mais 1h e então foram coletados 1 ml de cada tubo para realizar a extração de RNA (500 ul em cada tubo + 2 volumes de RNA protect).

Para o biofilme: 
-    O mesmo inóculo utilizado para as curvas de crescimento, foi utilizado para inocular uma placa de 12 poços. Após 24 h de incubação a 37 C, o crescimento (células planctônicas) foram transferidos para tubos para a leitura da OD600.
-    O biofilme aderido a placa foi suspenso em 3 ml de PBS (mesmo volume que o usado no início do experimento). Nessa suspensão também foi realizada a leitura da OD600.
-    Foram coletados 1m da suspensão do biofilme para a extração do RNA. Os RNAS foram extraídos com o kit: RNeasy Minikit da Qiagen.

Para os experimentos com as diferentes fontes de carbono: 
-    Foi utilizado o meio M9, e neste meio foram adicionados os açúcares a 40mM. A Pseudomonas foi inoculada e realizamos uma curva de crescimento. Conforme elas foram atingindo a mid-log, alíquotas foram sendo retiradas para extrair RNA.




  

