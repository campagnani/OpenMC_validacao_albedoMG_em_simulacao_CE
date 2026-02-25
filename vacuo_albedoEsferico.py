import openmc
import numpy as np
import os
from datetime import datetime


# Criar pasta e mudar para a pasta criada:
def mkdir(nome="teste_sem_nome",data=True,voltar=False):
    if (voltar==True):
        os.chdir("../")
    if (data==True):
        agora = datetime.now()
        nome = agora.strftime(nome+"_%Y%m%d_%H%M%S")
    if not os.path.exists(nome):
        os.makedirs(nome)
    os.chdir(nome)




# Mudar para a pasta
def chdir(nome=None):
    if (nome != None):
        os.chdir(nome)
    else:
        diretorio_atual = os.getcwd()
        diretorios = [diretorio for diretorio in os.listdir(diretorio_atual) if os.path.isdir(os.path.join(diretorio_atual, diretorio))]

        data_mais_recente = 0
        pasta_mais_recente = None

        for diretorio in diretorios:
            data_criacao = os.path.getctime(os.path.join(diretorio_atual, diretorio))
            if data_criacao > data_mais_recente:
                data_mais_recente = data_criacao
                pasta_mais_recente = diretorio

        if pasta_mais_recente:
            os.chdir(os.path.join(diretorio_atual, pasta_mais_recente))
            print("Diretório mais recente encontrado:", pasta_mais_recente)
        else:
            print("Não foi possível encontrar um diretório mais recente.")


mkdir(nome="simuAlbedoBasica", data=False, voltar=True)

# =============================================================================
# 1. DEFINIÇÃO DA MATRIZ DE ALBEDO (2 GRUPOS)
# =============================================================================
# Grupos: G1 (Rápido: 0.625 eV a 20 MeV), G2 (Térmico: 0.0 eV a 0.625 eV)
energia = [0.0, 0.625, 20.0e6]

# Matriz de probabilidade
# Linha 0 (Bate Térmico) -> 50% volta Térmico, 0% volta Rápido
# Linha 1 (Bate Rápido)  -> 80% volta Térmico, 10% volta Rápido
matriz = np.array([
    [0.5, 0.0], 
    [0.8, 0.1]  
])

vetor = np.array([
    [0.7, 0.2]  
])


# =============================================================================
# 2. GEOMETRIA COM A NOVA SUPERFÍCIE
# =============================================================================
# Usando a sua nova API!
esfera_refletora = openmc.Sphere(
    r=10.0, 
    boundary_type='reflective',
    albedo=matriz, 
    albedo_energy_grid=energia
)

# Célula preenchida com vácuo (para o nêutron viajar livre e bater na parede)
celula_interna = openmc.Cell(region=-esfera_refletora)
geometria = openmc.Geometry([celula_interna])

# =============================================================================
# 3. CONFIGURAÇÕES E FONTE
# =============================================================================
settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.batches = 100
settings.particles = 1000

# Fonte pontual no centro, emitindo nêutrons a 2 MeV (Grupo Rápido)
fonte = openmc.IndependentSource()
fonte.space = openmc.stats.Point((0.0, 0.0, 0.0))
fonte.energy = openmc.stats.Discrete([2.0e6], [1.0])
settings.source = fonte


# =============================================================================
# 4. TALLIES (O Espectrômetro)
# =============================================================================
tallies = openmc.Tallies()

# Cria um filtro de energia com os mesmos grupos da sua matriz
filtro_energia = openmc.EnergyFilter(energia)

# Cria um Tally para medir o fluxo dentro da nossa célula de vácuo
tally_fluxo = openmc.Tally(name='fluxo_interno')
tally_fluxo.filters = [openmc.CellFilter(celula_interna), filtro_energia]
tally_fluxo.scores = ['flux']
tallies.append(tally_fluxo)


# =============================================================================
# 5. EXPORTAÇÃO E EXECUÇÃO
# =============================================================================
geometria.export_to_xml()
settings.export_to_xml()
openmc.Materials().export_to_xml()
tallies.export_to_xml()

# Caminho absoluto ou relativo para o executável que você acabou de compilar
# Ajuste o caminho abaixo para onde o seu binário 'openmc' está!
caminho_executavel = os.path.abspath("../openmc/build/bin/openmc")

print(f"Rodando o OpenMC compilado localmente em: {caminho_executavel}")
openmc.run(openmc_exec=caminho_executavel)

# =============================================================================
# 6. PÓS-PROCESSAMENTO (Lendo a prova do crime)
# =============================================================================
print("\n" + "="*50)
print(" RESULTADOS DO ESPECTRO DENTRO DA ESFERA DE VÁCUO")
print("="*50)

# Abre o arquivo de estado gerado no final do batch 10
sp = openmc.StatePoint(f'statepoint.{settings.batches}.h5')

# Pega o nosso tally
tally = sp.get_tally(name='fluxo_interno')

# Imprime na tela usando o Pandas para ficar uma tabela bonita
df = tally.get_pandas_dataframe()

# Formata a impressão para vermos claramente os grupos
for index, row in df.iterrows():
    e_min = row['energy low [eV]']
    e_max = row['energy high [eV]']
    fluxo = row['mean']
    erro = row['std. dev.']
    
    grupo = "Térmico" if e_min == 0.0 else "Rápido"
    print(f"Grupo {grupo} ({e_min} a {e_max} eV): Fluxo = {fluxo:.5e} +/- {erro:.5e}")
