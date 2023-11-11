import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
# 環状化合物をきれいに描画 https://future-chem.com/rdkit-coordgen/
rdDepictor.SetPreferCoordGen(True)



# Header
st.title("Retro-inversian")
st.write("")

# query
# st.selectbox("Select the Constraint Type of the Peptide",
#              ['Liner'])
FASTA = st.text_area('input FASTA')
fasta = FASTA[::-1]

st.write('Retro-inversian replied "' + fasta.lower() + '"')

type = st.selectbox('Select the Constraint Type of the peptide', 
                    ['Liner', 'Head to Tail', 'Disulfide']
                    )

# %%
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdMolAlign, Draw, rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole
import streamlit as st
import py3Dmol
import copy

# aminoacid format
laa = 'N[C@@]([H])({X})C(=O)'
daa = 'N[C@@]({X})([H])C(=O)'

# residue dict
amino_acids = {
    'G': '[H]',
    'A': 'C',
    'V': 'C(C)C',
    'L': 'CC(C)C',
    'I': 'C(C)CC',
    'M': 'CCSC',
    'F': 'CC1=CC=CC=C1',
    'W': 'CC1=CNC2=C1C=CC=C2',
    'S': 'CO',
    'C': 'CS',
    'Y': 'CC1=CC=C(C=C1)O',
    'N': 'C(C(=O)N)',
    'Q': 'C(CC(=O)N)',
    'D': 'C(C(=O)O)',
    'E': 'C(CC(=O)O)',
    'H': 'C(C1=CNC=N1)',
    'K': 'CCCCN',
    'R': 'CCCNC(=N)N',
}

# insert variable to format
laaa = {k: laa.format(X=v) for k, v in amino_acids.items()}
daaa = {k: daa.format(X=v) for k, v in amino_acids.items()}

# add prolines and threonines
laaa['P'] = 'N1[C@@]([H])(CCC1)C(=O)'
daaa['P'] = 'N1[C@@](CCC1)([H])C(=O)'
laaa['T'] = 'N[C@@]([H])([C@]([H])(O)C)C(=O)'
daaa['T'] = 'N[C@@]([C@]([H])(O)C)([H])C(=O)'


FASTA = 'CAWAFAAAC'
fasta = FASTA[::-1]
type = 'Liner'



l_sec = ''.join([laaa[aa] for aa in FASTA])
d_sec = ''.join([daaa[aa] for aa in fasta])

# liner
l_liner = l_sec + 'O'
d_liner = d_sec + 'O'

l_ht = 'N9' + l_sec[1:-4] + '9(=O)'
d_ht = 'N9' + d_sec[1:-4] + '9(=O)'


# # Cyclization
# if '2' in l_sec:
#     l_ht = 'N3' + l_sec[1:-4] + '3(=O)'
#     laaa['J'] = 'N[C@@]([H])(CS3)C(=O)'
#     daaa['J'] = 'N[C@@](CS3)([H])C(=O)'
#     FASTA = FASTA.replace('C', 'J')
#     fasta = FASTA[::-1]
#     l_ds = ''.join([laaa[aa] for aa in FASTA])
#     d_ds = ''.join([daaa[aa] for aa in fasta])
# elif '1' in l_liner:
#     l_ht = 'N2' + l_sec[1:-4] + '2(=O)'
# else:
#     l_ht = 'N1' + l_sec[1:-4] + '1(=O)'
# if '2' in d_sec:
#     d_ht = 'N3' + d_sec[1:-4] + '3(=O)'
# elif '1' in l_liner:
#     d_ht = 'N2' + d_sec[1:-4] + '2(=O)'
# else:
#     d_ht = 'N1' + d_sec[1:-4] + '1(=O)'


# Draw.MolToImage(Chem.MolFromSmiles(l_ht), size=(1000, 400))

if type == 'Head to Tail':
    l_smiles = l_ht
    d_smiles = d_ht
elif type == 'Disulfide':
    l_smiles = l_ds
    d_smiles = d_ds
else:
    l_smiles = l_liner
    d_smiles = d_liner

# l_mol
l_mol = Chem.MolFromSmiles(l_smiles)
d_mol = Chem.MolFromSmiles(d_smiles)

# generate 3D conformations
molecule = Chem.AddHs(l_mol)
AllChem.EmbedMolecule(l_mol)
AllChem.UFFOptimizeMolecule(l_mol)
molecule = Chem.AddHs(d_mol)
AllChem.EmbedMolecule(d_mol)
AllChem.UFFOptimizeMolecule(d_mol)
mols = [l_mol, d_mol]



# %%
# output
st.image(Draw.MolToImage(l_mol, size=(1000, 400)))
st.code(Chem.MolToSmiles(l_mol))
st.image(Draw.MolToImage(d_mol, size=(1000, 400)))
st.code(Chem.MolToSmiles(d_mol))

# saves
col1, col2, col3 = st.columns(3)
with col1:
    bt1 = st.button('save L.sdf')
with col2:
    bt2 = st.button('save D.sdf')

if bt1:
    # 3D座標を生成
    molecule = Chem.AddHs(l_mol)
    AllChem.EmbedMolecule(l_mol)
    AllChem.UFFOptimizeMolecule(l_mol)

    # SDFファイルに保存
    writer = Chem.SDWriter(FASTA + '.sdf')
    writer.write(l_mol)
    writer.close()

if bt2:
    # 3D座標を生成
    molecule = Chem.AddHs(d_mol)
    AllChem.EmbedMolecule(d_mol)
    AllChem.UFFOptimizeMolecule(d_mol)

    # SDFファイルに保存
    writer = Chem.SDWriter(fasta.lower() + '.sdf')
    writer.write(d_mol)
    writer.close()


# %%

