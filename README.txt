Lumpac View

Programa que automatiza a construção de estruturas.


MODIFICACOES PROPOSTAS

--> A leitura do ligante acontece em uma etapa e sua inicializacao (x1 e x2)
    acontece em outra. Isso foi criado para facilitar o error handling.
    seria melhor fazer de uma vez só.

--> Quando ele não acha o arquivo no rmsd ele dá pau.
    precisa de um error handling.


Funcionamento

3 módulos

-> criacao de pontos.

-> ajuste do SA.

-> construcao dos complexos.



INICIO

1. Le o lumpacViewInput.txt e salva os dados no objeto readInput.
   - Procura os ligantes pelo nome, le as coordenadas e os salva para uso futuro.

2. Usa a matemática do Simas para definir os pontos 
   X1 - centro da coordenacao.
   X2 - direção da coordenacao.

3. ComplexCreator.start
   3.1. Ordena os ligantes: tridentados, bidentados e monodentados.
        !(seria necessario uma revisao em varios pontos para descrever polidentados)
   3.2. Obtencao dos pontos na esfera.
   3.3. Esticamento dos pontos.
   3.4. Gera posicoes iniciais usando os pontos e roda para o centro.

4. ComplexCreator.simulatedAnnealing
   -> Algoritmo padrao - a unica diferenca e uma temperatura variavel que inventei.
   -> calculateAllFit - passa a bola para o que esta implementado em Ligand.calculateFit.
   -> perturbOperations - 1. gera um vetor no centro e roda ao redor dele
                          2. roda em torno da molecula.
   -> No fim da simulacao retorna vector<CoordXYZ> allAtoms.

5. ControlMopac.optimize(AllAtoms)
   -> roda uma otimizacao do mopac com todos os atomos.
   -> se a otimizacao der certo, roda frequencia.

FIM



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
dentity NOME // monodentate, bidentate ou tridentate
label x y z  // coordendas no formato xyz
...
...
]
!!! IMPORTANTE
Monodentado - X1 = primeiro atomo.
              X2 = centro de geometrico -> X1.
              
Bidentado - X1 = (xyzAtomo1 + xyzAtomo2) / 2 
            X2 = Usa o terceiro atomo, varias operacoes.

Tridentado - X1 = baricentro dos 3 ligantes.
             X2 = normal do triangulo saindo da molecula.


Para parametrizar o simulated annealing ele lê
um ponto com esse formato:
point.txt
[
5       - variaveis
0.1e0   - temperature update
2.0     - max alpha angle
2.0     - max beta angle
500     - initial temperature
0.5     - acceptance probability
]

e entrega: 
fitness.txt
[
0.4214   - erro.
]


Rotinas:

AuxMath -> funções de álgebra linear.

Coordstructs -> struct usados no programa.

ControlMopac -> Controla a montagem dos inputs, execucao
                e analise dos resultados do mopac.

Fitness -> onde esta implementado a interacao entre os ligantes.

Ligand -> Cada objeto é um ligante todo.
          Centro, vetor direcional e ângulo com o lantanídeo.

MyExceptions -> Mensagens de erro.

Points -> onde estao os pontos da esfera.

ReadInput -> Le o LumpacViewInput.txt e gera os ligantes.
             Contem vector<Ligand> allLigands;

ReadMopac -> objeto usado pelo controlmopac para ler
             outputs do mopac.



WARNING

- em Ligand.cpp -> limite de atomos para bidentate e tridentate.
                   bi  precisa de 3 atomos
                   tri precisa de 4

- em ComplexCreator.perturbOperations -> ele mexe em todos os ligantes de uma vez so
                                         se ele mexe em um de cada vez poderia ser mais
                                         eficiente.

- em ControlMopac.buildMopacInput -> na frequencia parece que nao esta indo em precisao dupla.


- na contrucao do bidentado -> o terceiro atomo esta pra tras, isso pode nao
                               ser verdade em alguns casos.


- No kabsch algortithm -> ATENCAO - algoritmo completamente ineficiente.
                                    e necessario corrigir isso.



IDEIAS

Opcao new ligand - voce me da as coordenadas e a denticao.
                   eu monto um complexo cheio chelation = 8
                   com aquele ligante. rodo e recupero ele
                   deformado, que sera melhor para o uso
                   no lumpacview.



TRECHOS DE CODIGOS
:criacao de inputs
	/*
	int charge = 1;
	int coordination = 8;
	string ligandName = "DUCNAQ-OONO.xyz";
	string mopacExecPath = "M2009_Ln_Orbitals.exe";
	vector<string> options(5);
	options[0] = "mopac2009";
	options[1] = "DUCNAQ-OONO";
	options[2] = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n"; //freq = AUX THERMO FORCE
																		//adding charge
	string chargeString;
	stringstream convert;
	convert << charge;
	convert >> chargeString;
	chargeString = " CHARGE=" + chargeString + " ";
	// charge computed
	options[2] += " NOLOG GEO-OK SCFCRT=1.D-10" + chargeString;
	options[3] = "Eu_spk";
	options[4] = "Eu";
	BuildComplex bc_;
	bc_.makeComplexOptimizingInMopac(ligandName, coordination, charge, options, mopacExecPath);
	*/






