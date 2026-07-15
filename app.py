import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import tempfile
import os

st.set_page_config(page_title="🧬 Beast BioAnalyzer", page_icon="🧬", layout="wide")

st.title("🧬 BEAST MODE BIOINFORMATICS SUITE")
st.markdown("---")

# Sidebar Menu
st.sidebar.header("🔬 Tools")
tool = st.sidebar.selectbox("Choose Tool", [
    "DNA ↔ RNA Conversion",
    "RNA → Protein Translation",
    "Mutation Detection",
    "Codon Usage Analysis",
    "Sequence Alignment",
    "FASTA/FASTQ Parser"
])

st.sidebar.markdown("---")
st.sidebar.info("Made with ❤️ by a 12-year-old coder")

# ============= DNA ↔ RNA CONVERSION =============
if tool == "DNA ↔ RNA Conversion":
    st.header("🔄 DNA ↔ RNA Conversion")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("DNA → RNA")
        dna_input = st.text_area("Enter DNA Sequence:", placeholder="ATCGATCG...", key="dna")
        if dna_input:
            rna_output = dna_input.upper().replace('T', 'U')
            st.success(f"✅ RNA: {rna_output}")
            st.info(f"📏 Length: {len(rna_output)} bases")
    
    with col2:
        st.subheader("RNA → DNA")
        rna_input = st.text_area("Enter RNA Sequence:", placeholder="AUCGAUCG...", key="rna")
        if rna_input:
            dna_output = rna_input.upper().replace('U', 'T')
            st.success(f"✅ DNA: {dna_output}")
            st.info(f"📏 Length: {len(dna_output)} bases")

# ============= PROTEIN TRANSLATION =============
elif tool == "RNA → Protein Translation":
    st.header("🧬 RNA → Protein Translation")
    
    rna_seq = st.text_area("Enter RNA Sequence:", placeholder="AUGCGAUAA...", height=100)
    
    if rna_seq:
        try:
            seq_obj = Seq(rna_seq.upper())
            protein = seq_obj.translate(to_stop=False)
            
            st.success("✅ Protein Sequence Generated!")
            st.code(str(protein), language="text")
            
            st.info(f"📏 Amino Acids: {len(protein)}")
            
            st.subheader("📊 Amino Acid Composition")
            aa_counts = {}
            for aa in protein:
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
            
            for aa, count in sorted(aa_counts.items()):
                st.write(f"**{aa}**: {count}")
                
        except Exception as e:
            st.error(f"❌ Error: {str(e)}")

# ============= MUTATION DETECTION =============
elif tool == "Mutation Detection":
    st.header("🔬 Mutation Detection")
    
    col1, col2 = st.columns(2)
    
    with col1:
        original = st.text_area("Original Sequence:", placeholder="ATCGATCG...", key="orig")
    
    with col2:
        mutated = st.text_area("Mutated Sequence:", placeholder="ATCGTTCG...", key="mut")
    
    if original and mutated:
        st.subheader("🔍 Results")
        
        mutations = []
        for i in range(min(len(original), len(mutated))):
            if original[i].upper() != mutated[i].upper():
                orig_base = original[i].upper()
                mut_base = mutated[i].upper()
                
                purines = {'A', 'G'}
                pyrimidines = {'C', 'T', 'U'}
                
                if (orig_base in purines and mut_base in purines) or \
                   (orig_base in pyrimidines and mut_base in pyrimidines):
                    m_type = "Transition"
                else:
                    m_type = "Transversion"
                
                mutations.append(f"Position {i+1}: **{orig_base} → {mut_base}** ({m_type})")
        
        if len(original) != len(mutated):
            mutations.append(f"⚠️ Length Difference: {len(original)} vs {len(mutated)} bp")
        
        if mutations:
            st.warning(f"Found {len(mutations)} mutation(s):")
            for mut in mutations:
                st.markdown(f"- {mut}")
        else:
            st.success("✅ No mutations detected! Sequences are identical.")

# ============= CODON USAGE =============
elif tool == "Codon Usage Analysis":
    st.header("📊 Codon Usage Analysis")
    
    dna_seq = st.text_area("Enter DNA Sequence:", placeholder="ATGCGATCG...", height=150)
    
    if dna_seq:
        dna_seq = dna_seq.upper().replace('U', 'T')
        
        if len(dna_seq) % 3 != 0:
            st.warning(f"⚠️ Sequence length ({len(dna_seq)}) is not divisible by 3. Trailing bases will be ignored.")
        
        codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq) - (len(dna_seq) % 3), 3)]
        total = len(codons)
        
        codon_counts = {}
        for codon in codons:
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
        
        st.subheader(f"📈 Top Codons (Total: {total})")
        
        sorted_codons = sorted(codon_counts.items(), key=lambda x: x[1], reverse=True)
        
        for codon, count in sorted_codons[:10]:
            freq = (count / total) * 100
            st.write(f"**{codon}**: {count} ({freq:.2f}%)")

# ============= SEQUENCE ALIGNMENT =============
elif tool == "Sequence Alignment":
    st.header("🧬 Sequence Alignment")
    
    col1, col2 = st.columns(2)
    
    with col1:
        seq1 = st.text_area("Sequence 1:", placeholder="ATCG...", height=100, key="s1")
    
    with col2:
        seq2 = st.text_area("Sequence 2:", placeholder="ATCG...", height=100, key="s2")
    
    if seq1 and seq2:
        try:
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5
            
            alignments = aligner.align(seq1.upper(), seq2.upper())
            best = alignments[0]
            
            st.success(f"✅ Alignment Score: **{best.score}**")
            
            st.subheader("Aligned Sequences:")
            st.text(best.format())
            
        except Exception as e:
            st.error(f"❌ Error: {str(e)}")

# ============= FASTA/FASTQ PARSER =============
elif tool == "FASTA/FASTQ Parser":
    st.header("📄 File Parser")
    
    file_type = st.radio("Select File Type:", ["FASTA", "FASTQ"])
    
    uploaded_file = st.file_uploader(f"Upload {file_type} file", type=['fasta', 'fa', 'fastq', 'fq'])
    
    if uploaded_file:
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=f".{file_type.lower()}") as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                tmp_path = tmp_file.name
            
            if file_type == "FASTA":
                records = list(SeqIO.parse(tmp_path, "fasta"))
                st.success(f"✅ Found {len(records)} sequence(s)")
                
                for i, record in enumerate(records[:5]):
                    with st.expander(f"Sequence {i+1}: {record.id}"):
                        st.write(f"**Length:** {len(record.seq)} bp")
                        st.write(f"**Description:** {record.description}")
                        st.code(str(record.seq)[:200] + "...", language="text")
            
            else:
                records = list(SeqIO.parse(tmp_path, "fastq"))
                st.success(f"✅ Found {len(records)} read(s)")
                
                for i, record in enumerate(records[:5]):
                    quals = record.letter_annotations["phred_quality"]
                    avg_q = sum(quals) / len(quals)
                    
                    with st.expander(f"Read {i+1}: {record.id}"):
                        st.write(f"**Length:** {len(record.seq)} bp")
                        st.write(f"**Avg Quality:** {avg_q:.2f}")
                        st.write(f"**Status:** {'✅ Excellent (Q30+)' if avg_q >= 30 else '✅ Good (Q20+)' if avg_q >= 20 else '️ Poor'}")
            
            os.unlink(tmp_path)
            
        except Exception as e:
            st.error(f"❌ Error parsing file: {str(e)}")

# Footer
st.markdown("---")
st.markdown("### 🚀 About This Tool")
st.markdown("""
This bioinformatics suite is built for:
- DNA/RNA sequence analysis
- Protein translation
- Mutation detection
- Educational purposes

**Built with Streamlit + Biopython** 🧬
""")