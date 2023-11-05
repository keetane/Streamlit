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

st.write('Retro-inversian replied')
st.write(fasta.lower())

# %%
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import IPythonConsole
import streamlit as st

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
    'R': 'CCCNC(=N)N'
}

# insert variable to format
laaa = {k: laa.format(X=v) for k, v in amino_acids.items()}
daaa = {k: daa.format(X=v) for k, v in amino_acids.items()}

# add prolines and threonines
laaa['P'] = 'N1[C@@]([H])(CCC1)C(=O)'
daaa['P'] = 'N1([H])[C@@](CCC1)C(=O)'
laaa['T'] = 'N[C@@]([H])([C@]([H])(O)C)C(=O)'
daaa['T'] = 'N[C@@]([C@]([H])(O)C)([H])C(=O)'


# FASTA = 'FASTA'
# fasta = FASTA[::-1]

l_sec = ''.join([laaa[aa] for aa in FASTA])
l_liner = l_sec + 'O'
l_img = Chem.MolFromSmiles(l_liner)
d_sec = ''.join([daaa[aa] for aa in fasta])
d_liner = d_sec + 'O'
d_img = Chem.MolFromSmiles(d_liner)


# %%
# output
st.image(Draw.MolToImage(l_img, size=(1000, 400)))
st.image(Draw.MolToImage(d_img, size=(1000, 400)))


# saves
col1, col2 = st.columns(2)
with col1:
    bt1 = st.button('save L.sdf')
with col2:
    bt2 = st.button('save D.sdf')

if bt1:
    # 3D座標を生成
    molecule = Chem.AddHs(l_img)
    AllChem.EmbedMolecule(l_img)
    AllChem.UFFOptimizeMolecule(l_img)

    # SDFファイルに保存
    writer = Chem.SDWriter(FASTA + '.sdf')
    writer.write(l_img)
    writer.close()

if bt2:
    # 3D座標を生成
    molecule = Chem.AddHs(d_img)
    AllChem.EmbedMolecule(d_img)
    AllChem.UFFOptimizeMolecule(d_img)

    # SDFファイルに保存
    writer = Chem.SDWriter(fasta.lower() + '.sdf')
    writer.write(d_img)
    writer.close()

# %%
