import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import gc_content, molecular_weight
import tempfile
import os
import pandas as pd
import re
import numpy as np

st.set_page_config(page_title="Beast BioAnalyzer v5.0", page_icon="🧬", layout="wide")

st.title("🧬 BEAST MODE BIOINFORMATICS SUITE v5.0")
st.markdown("---")

# Sidebar Menu
st.sidebar.header(" Tools")
tool = st.sidebar.selectbox("Choose Tool", [
    "1. DNA ↔ RNA Conversion",
    "2. RNA → Protein Translation",
    "3. Mutation Detection",
    "4. Codon Usage Analysis",
    "5. Sequence Alignment (Global)",
    "6. FASTA/FASTQ Parser",
    "7. Reverse Complement",
    "8. GC Content & Melting Temp",
    "9. Restriction Enzyme Sites",
    "10. ORF Finder",
    "11. Nucleotide Frequency Chart",
    "12. Central Dogma (DNA->RNA->Protein)",
    "13. Hamming Distance",
    "14. Molecular Weight Calculator",
    "15. Motif Finder (Pattern Search)",
    "16. Advanced Graphing & Visualization",
    "17. CRISPR-Cas9 Cut Site & Efficiency",
    "18. Phylogenetic Distance Matrix",
    "19. PCR Primer Designer",
    "20. Translation Table Selector",
    "21. Sequence Statistics"
])

st.sidebar.markdown("---")
st.sidebar.info("Made with ❤️ by a 12-year-old coder")

# ============= 1. DNA ↔ RNA CONVERSION =============
if tool == "1. DNA ↔ RNA Conversion":
    st.header("🔄 DNA ↔ RNA Conversion")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("DNA → RNA")
        dna_input = st.text_area("Enter DNA Sequence:", placeholder="ATCGATCG...", key="dna1")
        if st.button("Convert to RNA", key="btn1"):
            if dna_input:
                rna_output = dna_input.upper().replace('T', 'U')
                st.success(f"RNA: {rna_output}")
                st.info(f"Length: {len(rna_output)} bases")
            else: st.warning("Please enter a sequence first.")
    with col2:
        st.subheader("RNA → DNA")
        rna_input = st.text_area("Enter RNA Sequence:", placeholder="AUCGAUCG...", key="rna1")
        if st.button("Convert to DNA", key="btn2"):
            if rna_input:
                dna_output = rna_input.upper().replace('U', 'T')
                st.success(f"DNA: {dna_output}")
                st.info(f"Length: {len(dna_output)} bases")
            else: st.warning("Please enter a sequence first.")

# ============= 2. PROTEIN TRANSLATION =============
elif tool == "2. RNA → Protein Translation":
    st.header("🧬 RNA → Protein Translation")
    rna_seq = st.text_area("Enter RNA Sequence:", placeholder="AUGCGAUAA...", height=100, key="rna2")
    table = st.selectbox("Genetic Code:", ["Standard", "Vertebrate Mitochondrial", "Invertebrate Mitochondrial", "Bacterial"])
    table_map = {"Standard": 1, "Vertebrate Mitochondrial": 2, "Invertebrate Mitochondrial": 5, "Bacterial": 11}
    
    if st.button("Translate to Protein", key="btn3"):
        if rna_seq:
            try:
                seq_obj = Seq(rna_seq.upper())
                protein = seq_obj.translate(table=table_map[table], to_stop=False)
                st.success("Protein Sequence Generated!")
                st.code(str(protein), language="text")
                st.info(f"Amino Acids: {len(protein)}")
            except Exception as e: st.error(f"Error: {str(e)}")
        else: st.warning("Please enter a sequence first.")

# ============= 3. MUTATION DETECTION =============
elif tool == "3. Mutation Detection":
    st.header("🔬 Mutation Detection")
    col1, col2 = st.columns(2)
    with col1: original = st.text_area("Original Sequence:", placeholder="ATCGATCG...", key="orig3")
    with col2: mutated = st.text_area("Mutated Sequence:", placeholder="ATCGTTCG...", key="mut3")
    if st.button("Detect Mutations", key="btn4"):
        if original and mutated:
            mutations = []
            for i in range(min(len(original), len(mutated))):
                if original[i].upper() != mutated[i].upper():
                    orig_base, mut_base = original[i].upper(), mutated[i].upper()
                    p, py = {'A', 'G'}, {'C', 'T', 'U'}
                    m_type = "Transition" if (orig_base in p and mut_base in p) or (orig_base in py and mut_base in py) else "Transversion"
                    mutations.append(f"Position {i+1}: {orig_base} → {mut_base} ({m_type})")
            if len(original) != len(mutated):
                mutations.append(f"Length Difference: {len(original)} vs {len(mutated)} bp")
            if mutations:
                st.warning(f"Found {len(mutations)} mutation(s):")
                for mut in mutations: st.markdown(f"- {mut}")
            else: st.success("No mutations detected!")
        else: st.warning("Please enter both sequences.")

# ============= 4. CODON USAGE =============
elif tool == "4. Codon Usage Analysis":
    st.header("📊 Codon Usage Analysis")
    dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATGCGATCG...", height=150, key="dna4")
    if st.button("Analyze Codons", key="btn5"):
        if dna_seq:
            dna_seq_clean = dna_seq.upper().replace('U', 'T')
            codons = [dna_seq_clean[i:i+3] for i in range(0, len(dna_seq_clean) - (len(dna_seq_clean) % 3), 3)]
            total = len(codons)
            codon_counts = {}
            for codon in codons: codon_counts[codon] = codon_counts.get(codon, 0) + 1
            
            st.subheader(f"Top Codons (Total: {total})")
            for codon, count in sorted(codon_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
                st.write(f"**{codon}**: {count} ({(count/total)*100:.2f}%)")
            
            df = pd.DataFrame(list(codon_counts.items()), columns=['Codon', 'Count'])
            df['Frequency(%)'] = (df['Count'] / total * 100).round(2)
            csv_data = df.to_csv(index=False).encode('utf-8')
            st.download_button(label="📥 Download Results as CSV", data=csv_data, file_name='codon_usage.csv', mime='text/csv', key='dl4')
        else: st.warning("Please enter a sequence first.")

# ============= 5. SEQUENCE ALIGNMENT =============
elif tool == "5. Sequence Alignment (Global)":
    st.header("🧬 Sequence Alignment")
    col1, col2 = st.columns(2)
    with col1: seq1 = st.text_area("Sequence 1:", placeholder="ATCG...", height=100, key="s1_5")
    with col2: seq2 = st.text_area("Sequence 2:", placeholder="ATCG...", height=100, key="s2_5")
    if st.button("Align Sequences", key="btn6"):
        if seq1 and seq2:
            try:
                aligner = PairwiseAligner()
                aligner.mode = 'global'; aligner.match_score = 2; aligner.mismatch_score = -1
                aligner.open_gap_score = -2; aligner.extend_gap_score = -0.5
                best = aligner.align(seq1.upper(), seq2.upper())[0]
                st.success(f"Alignment Score: **{best.score}**")
                st.text(best.format())
            except Exception as e: st.error(f"Error: {str(e)}")
        else: st.warning("Please enter both sequences.")

# ============= 6. FASTA/FASTQ PARSER =============
elif tool == "6. FASTA/FASTQ Parser":
    st.header("📄 File Parser")
    file_type = st.radio("Select File Type:", ["FASTA", "FASTQ"], key="radio6")
    uploaded_file = st.file_uploader(f"Upload {file_type} file", type=['fasta', 'fa', 'fastq', 'fq'], key="up6")
    if uploaded_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=f".{file_type.lower()}") as tmp:
            tmp.write(uploaded_file.getvalue()); tmp_path = tmp.name
        try:
            records = list(SeqIO.parse(tmp_path, file_type.lower()))
            st.success(f"Found {len(records)} sequence(s)")
            for i, record in enumerate(records[:5]):
                with st.expander(f"{file_type} {i+1}: {record.id}"):
                    st.write(f"Length: {len(record.seq)} bp")
                    if file_type == "FASTQ":
                        quals = record.letter_annotations["phred_quality"]
                        st.write(f"Avg Quality: {sum(quals)/len(quals):.2f}")
                        q30 = sum(1 for q in quals if q >= 30) / len(quals) * 100
                        st.write(f"Q30 Score: {q30:.2f}%")
        except Exception as e: st.error(f"Error: {str(e)}")
        os.unlink(tmp_path)

# ============= 7. REVERSE COMPLEMENT =============
elif tool == "7. Reverse Complement":
    st.header("🔄 Reverse Complement")
    dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATCG...", height=100, key="dna7")
    if st.button("Find Reverse Complement", key="btn7"):
        if dna_seq:
            try:
                seq_obj = Seq(dna_seq.upper())
                st.info(f"Original: {dna_seq.upper()}")
                st.success(f"Reverse Complement: {seq_obj.reverse_complement()}")
            except Exception as e: st.error(f"Error: {str(e)}")
        else: st.warning("Please enter a sequence first.")

# ============= 8. GC CONTENT & TM =============
elif tool == "8. GC Content & Melting Temp":
    st.header("️ GC Content & Melting Temperature")
    dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATCG...", height=100, key="dna8")
    if st.button("Calculate GC & Tm", key="btn8"):
        if dna_seq:
            dna_seq = dna_seq.upper().replace('U', 'T')
            length = len(dna_seq)
            gc_percentage = gc_content(dna_seq) * 100  # FIXED: Convert fraction to percentage
            gc_count = dna_seq.count('G') + dna_seq.count('C')
            tm = 64.9 + 41 * (gc_count - 16.4) / length if length >= 14 else (gc_count * 4) + ((length - gc_count) * 2)
            col1, col2, col3 = st.columns(3)
            with col1: st.metric("Length", f"{length} bp")
            with col2: st.metric("GC Content", f"{gc_percentage:.2f}%")
            with col3: st.metric("Melting Temp (Tm)", f"{tm:.1f}°C")
            if 40 <= gc_percentage <= 60: st.success("Optimal GC content (40-60%)")
            else: st.warning("Non-optimal GC content")
        else: st.warning("Please enter a sequence first.")

# ============= 9. RESTRICTION ENZYMES =============
elif tool == "9. Restriction Enzyme Sites":
    st.header("✂️ Restriction Enzyme Sites")
    dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATCG...", height=150, key="dna9")
    if st.button("Find Enzyme Sites", key="btn9"):
        if dna_seq:
            dna_seq = dna_seq.upper()
            enzymes = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "NotI": "GCGGCCGC", "XhoI": "CTCGAG", "PstI": "CTGCAG", "SalI": "GTCGAC"}
            found = False
            for enzyme, site in enzymes.items():
                if site in dna_seq:
                    positions = [i+1 for i in range(len(dna_seq)) if dna_seq.startswith(site, i)]
                    st.write(f"**{enzyme}** ({site}): Positions {positions}")
                    found = True
            if not found: st.info("No common restriction sites found")
        else: st.warning("Please enter a sequence first.")

# ============= 10. ORF FINDER =============
elif tool == "10. ORF Finder":
    st.header("🔍 Open Reading Frame (ORF) Finder")
    dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATGCGATCG...", height=150, key="dna10")
    min_orf = st.slider("Minimum ORF Length (aa):", 10, 100, 30)
    if st.button("Find ORFs", key="btn10"):
        if dna_seq:
            try:
                seq_obj = Seq(dna_seq.upper())
                orfs = seq_obj.translate(to_stop=False).split('*')
                valid_orfs = [orf for orf in orfs if len(orf) >= min_orf]
                st.write(f"Found {len(valid_orfs)} ORF(s) >= {min_orf} amino acids")
                for i, orf in enumerate(valid_orfs[:5], 1):
                    st.write(f"ORF #{i}: {len(orf)} amino acids")
                    st.code(str(orf)[:50] + "...", language="text")
            except Exception as e: st.error(f"Error: {str(e)}")
        else: st.warning("Please enter a sequence first.")

# ============= 11. NUCLEOTIDE FREQUENCY =============
elif tool == "11. Nucleotide Frequency Chart":
    st.header("📊 Nucleotide Frequency Chart")
    dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATCG...", height=100, key="dna11")
    if st.button("Generate Chart", key="btn11"):
        if dna_seq:
            dna_seq = dna_seq.upper()
            counts = {'A': dna_seq.count('A'), 'T': dna_seq.count('T'), 'G': dna_seq.count('G'), 'C': dna_seq.count('C')}
            df = pd.DataFrame(list(counts.items()), columns=['Base', 'Count']).set_index('Base')
            st.bar_chart(df)
        else: st.warning("Please enter a sequence first.")

# ============= 12. CENTRAL DOGMA =============
elif tool == "12. Central Dogma (DNA->RNA->Protein)":
    st.header("🧬 Central Dogma (DNA → RNA → Protein)")
    dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATGCGATCG...", height=100, key="dna12")
    if st.button("Run Central Dogma", key="btn12"):
        if dna_seq:
            dna_seq = dna_seq.upper()
            rna = dna_seq.replace('T', 'U')
            protein = Seq(rna).translate(to_stop=False)
            st.info(f"**DNA:** {dna_seq}")
            st.success(f"**RNA:** {rna}")
            st.write(f"**Protein:** {protein}")
        else: st.warning("Please enter a sequence first.")

# ============= 13. HAMMING DISTANCE =============
elif tool == "13. Hamming Distance":
    st.header("📏 Hamming Distance Calculator")
    col1, col2 = st.columns(2)
    with col1: s1 = st.text_area("Sequence 1:", placeholder="ATCG", height=100, key="h1")
    with col2: s2 = st.text_area("Sequence 2:", placeholder="ATCG", height=100, key="h2")
    if st.button("Calculate Distance", key="btn13"):
        if s1 and s2:
            if len(s1) != len(s2): st.error("Sequences must be of equal length!")
            else:
                distance = sum(c1 != c2 for c1, c2 in zip(s1.upper(), s2.upper()))
                st.success(f"Hamming Distance: **{distance}**")
                st.info(f"Similarity: {100 - (distance / len(s1) * 100):.2f}%")
        else: st.warning("Please enter both sequences.")

# ============= 14. MOLECULAR WEIGHT =============
elif tool == "14. Molecular Weight Calculator":
    st.header("⚖️ Molecular Weight Calculator")
    seq_input = st.text_area("Enter Sequence:", placeholder="ATCG or Protein...", height=100, key="mw1")
    seq_type = st.selectbox("Select Type:", ["DNA", "RNA", "Protein"], key="mw_type")
    if st.button("Calculate Weight", key="btn14"):
        if seq_input:
            length = len(seq_input.strip())
            weights = {"DNA": 650, "RNA": 340, "Protein": 110}
            st.success(f"Approx. Molecular Weight: **{length * weights[seq_type]:,} Daltons (Da)**")
        else: st.warning("Please enter a sequence first.")

# ============= 15. MOTIF FINDER =============
elif tool == "15. Motif Finder (Pattern Search)":
    st.header(" Motif Finder")
    main_seq = st.text_area("Main Sequence:", placeholder="ATCGATCG...", height=100, key="main_motif")
    motif = st.text_input("Motif to find:", placeholder="ATCG", key="motif_search")
    if st.button("Find Motif", key="btn15"):
        if main_seq and motif:
            main_seq, motif = main_seq.upper(), motif.upper()
            positions = [i+1 for i in range(len(main_seq)) if main_seq.startswith(motif, i)]
            if positions: st.success(f"Found {len(positions)} occurrence(s) at positions: {positions}")
            else: st.info("Motif not found.")
        else: st.warning("Please enter both sequences.")

# ============= 16. ADVANCED GRAPHING =============
elif tool == "16. Advanced Graphing & Visualization":
    st.header(" Advanced Graphing & Visualization")
    graph_type = st.selectbox("Select Graph Type:", ["Amino Acid Composition", "GC Content Sliding Window", "Nucleotide Distribution"])
    
    if graph_type == "Amino Acid Composition":
        protein_seq = st.text_area("Enter Protein Sequence:", placeholder="MKTIIALSY...", height=100, key="aa_graph")
        if st.button("Generate Composition Chart", key="btn_aa"):
            if protein_seq:
                protein_seq = protein_seq.upper().replace('*', '')
                aa_counts = {}
                for aa in protein_seq: aa_counts[aa] = aa_counts.get(aa, 0) + 1
                df = pd.DataFrame(list(aa_counts.items()), columns=['Amino_Acid', 'Count']).sort_values('Count', ascending=False).set_index('Amino_Acid')
                st.bar_chart(df)
    
    elif graph_type == "GC Content Sliding Window":
        dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATCGATCG...", height=150, key="gc_graph")
        window_size = st.slider("Window Size:", 10, 100, 50)
        if st.button("Generate GC Graph", key="btn_gc"):
            if dna_seq:
                dna_seq = dna_seq.upper().replace('U', 'T')
                gc_values, positions = [], []
                for i in range(0, len(dna_seq) - window_size + 1, 10):
                    window = dna_seq[i:i+window_size]
                    gc_percentage = gc_content(window) * 100  # FIXED: Convert fraction to percentage
                    gc_values.append(gc_percentage)
                    positions.append(i + window_size//2)
                st.line_chart(pd.DataFrame({'GC_Content': gc_values}, index=positions))
    
    elif graph_type == "Nucleotide Distribution":
        dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATCG...", height=100, key="nucl_graph")
        if st.button("Generate Distribution", key="btn_nucl"):
            if dna_seq:
                seq = dna_seq.upper()
                counts = {'A': seq.count('A'), 'T': seq.count('T'), 'G': seq.count('G'), 'C': seq.count('C')}
                df = pd.DataFrame(list(counts.items()), columns=['Base', 'Count']).set_index('Base')
                st.bar_chart(df)

# ============= 17. CRISPR-CAS9 PREDICTOR =============
elif tool == "17. CRISPR-Cas9 Cut Site & Efficiency":
    st.header(" CRISPR-Cas9 Cut Site & Efficiency Predictor")
    st.markdown("Predicts SpCas9 cut sites (3bp upstream of NGG PAM) and calculates cleavage efficiency.")
    
    target_dna = st.text_area("Enter Target DNA Sequence:", placeholder="ATCGATCG...", height=150, key="crispr_dna")
    
    if st.button("🔍 Predict CRISPR Cut Sites", key="btn_crispr"):
        if target_dna:
            seq = target_dna.upper().replace('U', 'T')
            regex_pam = ".GG"
            matches = [(m.start(), m.group()) for m in re.finditer(f'(?={regex_pam})', seq)]
            
            if not matches:
                st.warning("No NGG PAM sites found in the sequence.")
            else:
                st.success(f"Found {len(matches)} potential PAM site(s).")
                
                for idx, (pos, pam) in enumerate(matches):
                    if pos >= 20:
                        protospacer = seq[pos-20:pos]
                        score = 50
                        gc_count = protospacer.count('G') + protospacer.count('C')
                        gc_pct = (gc_count / 20) * 100
                        
                        if 40 <= gc_pct <= 60: score += 20
                        elif 30 <= gc_pct <= 70: score += 10
                        if protospacer[19] in ['G', 'C']: score += 15
                        if protospacer[0] == 'G': score -= 10
                        if 'TTTT' in protospacer: score -= 25
                        score = max(0, min(100, score))
                        
                        tier = "HIGH" if score >= 70 else "MODERATE" if score >= 40 else "LOW"
                        cut_pos = pos - 3
                        
                        with st.expander(f"Site #{idx+1}: Position {pos+1} (Score: {score}% - {tier})"):
                            st.markdown(f"**Protospacer:** `{protospacer}`")
                            st.markdown(f"**PAM:** `{pam}`")
                            st.markdown(f"**Cut Site:** Between bp {cut_pos} and {cut_pos+1}")
                            st.write(f"**GC Co