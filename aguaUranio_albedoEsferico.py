import openmc
import numpy as np
import os
from datetime import datetime
import sys

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

ROOT_DIR = os.getcwd()

CAMINHO_EXEC = os.path.abspath("/home/thalles/0-git/openmc2/build0/bin/openmc")

arquivo_albedo = os.path.join(ROOT_DIR, "albedos.py")

if os.path.exists(arquivo_albedo):
    print("="*60)
    print(" ARQUIVO 'albedos.py' ENCONTRADO! IMPORTANDO DADOS...")
    print("="*60)
    # Adiciona o diretório atual ao path do Python para importar o arquivo
    sys.path.append(ROOT_DIR)
    import albedos
    
    matriz_albedo = albedos.matriz_albedo
    limites_energia = albedos.limites_energia
    num_grupos = len(limites_energia) - 1

else:
    print("="*60)
    print(" ARQUIVO 'albedos.py' NÃO ENCONTRADO. INICIANDO GERAÇÃO...")
    print("="*60)
    mkdir(nome="genAlbedo", data=False, voltar=True)

    # 1. MATERIAIS E GEOMETRIA DO REFLETOR
    agua = openmc.Material(name='Agua Leve')
    agua.add_nuclide('H1', 2.0)
    agua.add_nuclide('O16', 1.0)
    #agua.add_s_alpha_beta('c_H_in_H2O')
    agua.set_density('g/cm3', 1.0)

    materiais = openmc.Materials([agua])
    materiais.export_to_xml()

    s_vacuo_interna = openmc.Sphere(r=29.0, boundary_type='vacuum')
    s_albedo = openmc.Sphere(r=30.0)
    s_vacuo_externa = openmc.Sphere(r=1000.0, boundary_type='vacuum')

    c_fonte = openmc.Cell(name='Vacuo fonte', fill=None, region=+s_vacuo_interna & -s_albedo)
    c_agua = openmc.Cell(name='Agua Refletora', fill=agua, region=+s_albedo & -s_vacuo_externa)

    geometria = openmc.Geometry([c_fonte, c_agua])

    # 2. DEFINIÇÃO DOS GRUPOS
    limites_energia = [0.0, 0.625, 5.0, 1.0e3, 100.0e3, 1.0e6, 20.0e6]
    num_grupos = len(limites_energia) - 1
    matriz_albedo = np.zeros((num_grupos, num_grupos))

    # 3. LOOP DE SIMULAÇÕES
    for g_in in range(num_grupos):
        E_min = limites_energia[g_in]
        E_max = limites_energia[g_in + 1]
        
        print(f"Rodando Simulação {g_in+1}/{num_grupos} -> Injetando nêutrons no Grupo {g_in} ({E_min} a {E_max} eV)...")
        
        settings = openmc.Settings()
        settings.run_mode = 'fixed source'
        settings.particles = 50000 
        settings.batches = 10
        
        fonte = openmc.IndependentSource()
        fonte.space = openmc.stats.SphericalIndependent(
            r=openmc.stats.Discrete([29.5], [1.0]),
            cos_theta=openmc.stats.Uniform(0.0, np.pi),
            phi=openmc.stats.Uniform(0.0, 2.0 * np.pi)
        )
        fonte.energy = openmc.stats.Uniform(E_min, E_max)
        settings.source = fonte
        
        tallies = openmc.Tallies()
        filtro_energia = openmc.EnergyFilter(limites_energia)
        filtro_superficie = openmc.SurfaceFilter([s_albedo])
        
        tally_ida = openmc.Tally(name='corrente_ida')
        tally_ida.filters = [filtro_superficie, openmc.CellFromFilter([c_fonte]), filtro_energia]
        tally_ida.scores = ['current']
        tallies.append(tally_ida)

        tally_volta = openmc.Tally(name='corrente_volta')
        tally_volta.filters = [filtro_superficie, openmc.CellFromFilter([c_agua]), filtro_energia]
        tally_volta.scores = ['current']
        tallies.append(tally_volta)
        
        materiais.export_to_xml()
        geometria.export_to_xml()
        settings.export_to_xml()
        tallies.export_to_xml()
        
        openmc.run(openmc_exec=CAMINHO_EXEC, output=False)
        
        # Extração de Dados
        sp = openmc.StatePoint(f'statepoint.{settings.batches}.h5')
        df_ida = sp.get_tally(name='corrente_ida').get_pandas_dataframe()
        df_volta = sp.get_tally(name='corrente_volta').get_pandas_dataframe()
        
        J_ida = 0.0
        for idx, row in df_ida.iterrows():
            if row['energy low [eV]'] == E_min:
                J_ida = abs(row['mean']) 
                break
                
        for idx, row in df_volta.iterrows():
            e_low = row['energy low [eV]']
            try:
                g_out = limites_energia.index(e_low)
            except ValueError:
                continue
                
            if g_out < num_grupos:
                J_volta = abs(row['mean'])
                probabilidade = (J_volta / J_ida) if J_ida > 0 else 0.0
                matriz_albedo[g_in, g_out] = probabilidade

        sp.close()
        os.remove(f'statepoint.{settings.batches}.h5') 

    print("\nMATRIZ DE ALBEDO OBTIDA COM SUCESSO!")
    
    # SALVANDO O ARQUIVO albedos.py NA RAIZ
    with open(arquivo_albedo, "w") as f:
        f.write("import numpy as np\n\n")
        f.write(f"limites_energia = {limites_energia}\n\n")
        f.write("matriz_albedo = np.array([\n")
        for linha in matriz_albedo:
            f.write("    [" + ", ".join([f"{val:.5f}" for val in linha]) + "],\n")
        f.write("])\n")
    print(f"Matriz salva no arquivo: {arquivo_albedo}\n")













# =============================================================================
# RE-DECLARAÇÃO DOS MATERIAIS PARA AS PRÓXIMAS ETAPAS
# =============================================================================
agua = openmc.Material(name='Agua Leve')
agua.add_nuclide('H1', 2.0)
agua.add_nuclide('O16', 1.0)
#agua.add_s_alpha_beta('c_H_in_H2O')
agua.set_density('g/cm3', 1.0)

uranio = openmc.Material(name='Uranio')
uranio.add_nuclide("U235", 1)
uranio.set_density('g/cm3', 18.0)

materiais = openmc.Materials([agua, uranio])

# =============================================================================
# SIMULAÇÃO COMPLETA (Padrão de Referência)
# =============================================================================
mkdir(nome="simuCompleta", data=False, voltar=True)

print("="*60)
print(" Simulação Completa")
print("="*60)

s_uranio = openmc.Sphere(r=15.0)
s_albedo_referencia = openmc.Sphere(r=30.0)
s_vacuo_externa_ref = openmc.Sphere(r=1000.0, boundary_type='vacuum')

c_uranio = openmc.Cell(name='Uranio', fill=uranio, region=-s_uranio)
c_agua_buffer = openmc.Cell(name='Agua Buffer', fill=agua, region=+s_uranio & -s_albedo_referencia)
c_agua = openmc.Cell(name='Agua Refletora', fill=agua, region=+s_albedo_referencia & -s_vacuo_externa_ref)

geometria = openmc.Geometry([c_uranio, c_agua_buffer, c_agua])

settings = openmc.Settings()
settings.particles = 10000
settings.batches = 80
settings.inactive = 10
settings.source = openmc.IndependentSource(space=openmc.stats.Point())
settings.output = {'tallies': False}

materiais.export_to_xml()
geometria.export_to_xml()
settings.export_to_xml()

openmc.run(openmc_exec=CAMINHO_EXEC)

# =============================================================================
# SIMULAÇÃO COM MATRIZ DE ALBEDO
# =============================================================================
mkdir(nome="simuAlbedo", data=False, voltar=True)

print("\n" + "="*60)
print(" Simulação Albedo (Com Matriz de " + str(num_grupos) + " grupos)")
print("="*60)

# A MÁGICA ACONTECE AQUI
s_albedo_real = openmc.Sphere(r=30.0, boundary_type='reflective', albedo=matriz_albedo, albedo_energy_grid=limites_energia)
s_uranio_albedo = openmc.Sphere(r=15.0)

c_uranio_albedo = openmc.Cell(name='Uranio', fill=uranio, region=-s_uranio_albedo)
c_agua_buffer_albedo = openmc.Cell(name='Agua Buffer', fill=agua, region=+s_uranio_albedo & -s_albedo_real)

geometria = openmc.Geometry([c_uranio_albedo, c_agua_buffer_albedo])

settings = openmc.Settings()
settings.particles = 10000 # Igualado a completa para comparação justa de tempo
settings.batches = 80
settings.inactive = 10
settings.source = openmc.IndependentSource(space=openmc.stats.Point())
settings.output = {'tallies': False}

materiais.export_to_xml()
geometria.export_to_xml()
settings.export_to_xml()

openmc.run(openmc_exec=CAMINHO_EXEC)