import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import GC, molecular_weight
import pandas as pd
import re

st.set_page_config(page_title="Titan BioAnalyzer", page_icon="", layout="wide")
st.title(" TITAN BIOINFORMATICS SUITE")
st.sidebar.header("🔬 Tools")
tool = st.sidebar.selectbox("Choose Tool", [
    "1. DNA ↔ RNA", "2. Translation", "3. Reverse Complement",
    "4. GC Content", "5. Restriction Sites", "6. Hamming Distance", "7. Motif Finder"
])
st.sidebar.info("© 2026 Titan Bioinformatics | Shivay Singh")

if tool == "1. DNA ↔ RNA":
    st.header(" DNA ↔ RNA Conversion")
    col1, col2 = st.columns(2)
    with col1:
        dna = st.text_area("DNA:", placeholder="ATCG...", key="d1")
        if st.button("DNA→RNA", key="b1"):
            if dna: st.success(f"RNA: {dna.upper().replace('T','U')}")
    with col2:
        rna = st.text_area("RNA:", placeholder="AUCG...", key="r1")
        if st.button("RNA→DNA", key="b2"):
            if rna: st.success(f"DNA: {rna.upper().replace('U','T')}")

elif tool == "2. Translation":
    st.header("🧬 RNA → Protein")
    rna = st.text_area("RNA:", placeholder="AUGCGA...", height=100, key="t1")
    if st.button("Translate", key="b3"):
        if rna:
            protein = Seq(rna.upper()).translate(to_stop=False)
            st.success(f"Protein: {protein}")
            st.info(f"Amino Acids: {len(protein)}")

elif tool == "3. Reverse Complement":
    st.header("🔄 Reverse Complement")
    dna = st.text_area("DNA:", placeholder="ATCG...", height=100, key="rc1")
    if st.button("Find", key="b4"):
        if dna: st.success(f"Result: {Seq(dna.upper()).reverse_complement()}")

elif tool == "4. GC Content":
    st.header("🌡️ GC Content & Melting Temp")
    dna = st.text_area("DNA:", placeholder="ATCG...", height=100, key="gc1")
    if st.button("Calculate", key="b5"):
        if dna:
            dna = dna.upper().replace('U','T')
            gc = GC(dna)
            gc_count = dna.count('G') + dna.count('C')
            length = len(dna)
            tm = 64.9 + 41*(gc_count - 16.4)/length if length >= 14 else gc_count*4 + (length-gc_count)*2
            col1, col2, col3 = st.columns(3)
            with col1: st.metric("Length", f"{length} bp")
            with col2: st.metric("GC %", f"{gc:.2f}%")
            with col3: st.metric("Tm", f"{tm:.1f}°C")

elif tool == "5. Restriction Sites":
    st.header("✂️ Restriction Enzyme Sites")
    dna = st.text_area("DNA:", placeholder="ATCG...", height=150, key="re1")
    if st.button("Find Sites", key="b6"):
        if dna:
            dna = dna.upper()
            enzymes = {"EcoRI":"GAATTC","BamHI":"GGATCC","HindIII":"AAGCTT","NotI":"GCGGCCGC"}
            for enz, site in enzymes.items():
                if site in dna:
                    pos = [i+1 for i in range(len(dna)) if dna.startswith(site, i)]
                    st.write(f"**{enz}** ({site}): Positions {pos}")

elif tool == "6. Hamming Distance":
    st.header("📏 Hamming Distance")
    col1, col2 = st.columns(2)
    with col1: s1 = st.text_area("Seq 1:", placeholder="ATCG", height=100, key="h1")
    with col2: s2 = st.text_area("Seq 2:", placeholder="ATCG", height=100, key="h2")
    if st.button("Calculate", key="b7"):
        if s1 and s2:
            if len(s1) != len(s2): st.error("Equal length required!")
            else:
                d = sum(a!=b for a,b in zip(s1.upper(),s2.upper()))
                st.success(f"Distance: {d} | Similarity: {100-(d/len(s1)*100):.2f}%")

elif tool == "7. Motif Finder":
    st.header(" Motif Finder")
    seq = st.text_area("Sequence:", placeholder="ATCGATCG...", height=100, key="m1")
    motif = st.text_input("Motif:", placeholder="ATCG", key="m2")
    if st.button("Find", key="b8"):
        if seq and motif:
            seq, motif = seq.upper(), motif.upper()
            pos = [i+1 for i in range(len(seq)) if seq.startswith(motif, i)]
            if pos: st.success(f"Found {len(pos)} times at: {pos}")
            else: st.info("Not found.")

st.markdown("---")
st.markdown("### 🚀 Titan Bioinformatics | © 2026 Shivay Singh")