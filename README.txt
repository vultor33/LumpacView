Lumpac View

Programa que automatiza a construção de estruturas.

Funcionamento

1. Le o lumpacViewInput.txt e salva os dados no objeto readInput.
   - Procura os ligantes pelo nome, le as coordenadas e os salva para uso futuro.

2. Usa a matemática do Simas para definir os pontos 
   X1 - centro da coordenacao.
   X2 - direção da coordenacao -> pX1 -pX2.

Tarefas

- Encontrar uma forma de tirar os ligantes do readInput.

FORMATO DOS INPUTS

LumpacViewInput.txt
[
Eu             // metal 
spk_Gd.inp     // nome do arquivo de parametros
6              // numero total de ligantes
h2o            // nome dos arquivos dos ligantes sem .xyz
btfa
tta
...
...
]

Ligante.xyz
[
n            // numero de atomos
dentity      // monodentate, bidentate ou tridentate
label x y z  // coordendas no formato xyz
...
...
]
!!! IMPORTANTE
Monodentado - ligacao definida com o primeiro atomo e os outros
              usados na direcao.
Bidentado - ligacao definida com os dois primeiros e no
            plano do terceiro.
Tridentado - ligacao definida com os tres primeiros. Os outros
             foram usados na direcao.


Rotinas:

AuxMath -> funções de álgebra linear.

Coordstructs -> struct usados no programa.

Ligand -> Cada objeto é um ligante todo.
          Centro, vetor direcional e ângulo com o lantanídeo.

MyExceptions -> Mensagens de erro.

ReadInput -> Le o LumpacViewInput.txt e gera os ligantes.
             Contem vector<Ligand> allLigands;





WARNING

- em Ligand.cpp -> limite de atomos para bidentate e tridentate.
