#################################################################################################################
#                                                 Parâmetros                                                    #
#################################################################################################################

N_max=5            # Número de iterações
E_corte=5         # Energia de Corte
G=1.0              # Lambda
Ed=-0.0           # Energia do nível unicamente ocupado
U=0.0             # Energia de respulsão
V=0.1         # Gamma   (proporcional ao quadrado da hibridização)

#################################################################################################################
#                                                 Importing                                                      #
#################################################################################################################

using LinearAlgebra  #Pacote de Algebra linear
using SpecialFunctions  #Pacote de Funções Especiais
using BitBasis, Base.Threads      # numbers in binary representation