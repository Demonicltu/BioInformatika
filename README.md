# BioInformatika

Laboratorinio darbo analizė bei funkcijų įgyvendintų kode aprašymas laisva forma.

# Trečias laboratorinis darbas

Pastabos: prieš leidžiant projektą reikia pakeisti 

BioInformatika/3 Laboratorinis darbas/NEW.idea aplanką į .idea

Dabar taip yra padaryta, nes kertasi trečiojo laboratorinio darbo .idea aplankas su pirmo laboratorinio darbo .idea

1) Apibūdinkite fastq formatą. (https://en.wikipedia.org/wiki/FASTQ_format). Kokia papildoma informacija pateikiam lyginant su FASTA formatu?

fasta ir fastQ formatų panašumai:
  - Abu formatai kaupia sekos duomenis bei sekos metaduomenis (angl. data and metadata).
  - Abu formatai yra pagrįsti tekstu (angl. text-based).

fasta ir fastQ formatų skirtumai:
 - fastQ formatas buvo sugalvotas išspręsti konkrečią sekos metu kylančią problemą: atsižvelgiant į tai, kaip veikia skirtingos sekvenavimo technologijos, tikimybė, kad teisingai identifikuotas bus nukleotidas labai skiriasi. 
 - fasta formatas neturi standartizuoto būdo tai užkuduoti. Užtat fastQ yra papildomai atsirandanti nukleotido kokybės balų seka. (angl. sequence of quality scores)
 - Paprasčiau šnekant, fastQ formatas papildomai padeda nustatyti sekos kokybę.
 - fastQ formato papildoma informacija: 

    a) Eilutė, prasidedanti @, kurioje yra sekos ID.

    b) Viena ar daugiau eilučių, kuriose yra seka.

    c) Nauja eilutė, prasidedanti simboliu +, tuščia arba kartojanti sekos ID.

    d) Viena ar daugiau eilučių, kuriose yra kokybės balai.

--------------------------------------------------------------------------------------------------------------------

2) Kurią mėnesio dieną Jūs gimėte? Prie dienos pridėkite 33. Koks ASCII simbolis atitinka šį skaičių?

Gimiau 30 dieną. Prie 30 pridėjus 33 gaunu 63, ASCII simbolio atitiktis (tarp kabučių ") yra "?"

--------------------------------------------------------------------------------------------------------------------

3) Kodėl pirmi 32 ASCII kodai negali būti naudojami sekos kokybei koduoti?

Jei būtų naudojami pirmi 32 ASCII simboliai, tai sekos galėtų įgauti neigiamą kokybės balą.

--------------------------------------------------------------------------------------------------------------------

4. Parašykite skriptą, kuris:
  a) nustatytu koks kokybės kodavimas yra naudojamas pateiktame faile. Galimos koduotės:
  i. Sanger Phred+33
  ii. Solexa Solexa+64
  iii. Illumina 1.3+ Phred+64
  iv. Illumina 1.5+ Phred+64
  v. Illumina 1.8+ Phred+33
  Parašykite, kokią koduotę nustatėte ir kuo remiantis?

Nustatytos Sanger Phread+33 arba Illumina 1.8+ Phread+33 koduotės, galutinai negalima pasakyti tikslios koduotės, nes Sanger Phread+33 ir Illumina 1.8+ Phread+3 yra labai panašios, jos skiriasi tik vienu ASCII simboliu "j", nuordytame faile nebuvo tokio simbolio, todėl negalima tiksliai nustatyti kokia iš šių dviejų koduočių buvo naudojama. Kitos koduotės netiko, nes buvo faile naudojami ASCII simboliai, kurie neįeina į atitinkamas koduotes.

--------------------------------------------------------------------------------------------------------------------

4.b) Analizuotų C/G nukleotidų pasiskirstymą read’uose. Pateikite grafiką, kurio y ašyje būtų read’ų skaičius, x ašyje - C/G nukletidų dalis read’o sekoje (100 proc. Reikštų, kad visi simboliai read’o sekoje yra G ir C) Parašykite, koks „stambių“ pikų skaičius yra gautame grafike? (tikrai mažiau nei 6)

Atlikau papildomus skaičiavimus, kad surasti C/G nukleotidų pasiskirstymą read'e. 
Buvo naudojama funkcija atlikti papildomiems skaičiavimams beieškant C/G nukleotidų pasiskirstymo, tai būtų: GCratio(sequence), kuri priima DNR seką ir grąžina sekoje esančių C/G nukleotidų santykį. Skaičiai apvalinami iki dviejų skaičių po kablelio tikslumo. 

Gauti rezultatai surašyti į excel dokumentą /3 Laboratorinis darbas/CG nukleotidų pasiskirstymo analizė.xlsx. Sugeneruotoje diagramoje matomi keli peak'ai, tai yra nuo 0.24 iki 0.39, nuo 0.47 iki 0.58 bei nuo 0.64 iki 0.76.

--------------------------------------------------------------------------------------------------------------------

4.c) Paimtų po 5 kiekvieno piko viršūnės sekų ir atliktų blast’o paieškas. Naudokite nr/nt duombazę, paiešką apribokite taip, kad ieškotų atitikmenų tik bakterinės sekose (organizmas “bacteria”). Analizei naudokite tik patį pirmą atitikmenį. Pateikite lentelę, kurioje būtų read’o id ir rasto mikroorganizmo rūšis

Buvo naudojama duombazė (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&BLAST_SPEC=MicrobialGenomes) atrasti bakterijas esančias mėginyje.

Pirmas peak'as:

Sekos GAAGTTCTAATGATTTATCGACTGCTAGACAGCAGATATGCTTAATTTAAGTACCTTATGCTAATTGGCACCATGGGAGTGGGACAGAAATGATATTTTCGTAAAATTTATTTCGTAGTCCCACCACAAATCGCATTGCCTGTAGAATTTC rezultatas => Stafilokokas – patogeninis mikroorganizmas. 
Detaliau galima rasti /3 Laboratorinis darbas/Peakai/Peakas 1 dalis 1.txt

Sekos AACTTCACTCTTTGGTATATTGTTATTAGCAGGAAAACGATAAAACCTTTAATCAGAGATATTGAACTTCGTTTCTTTCTGTTAATAGCCTTAGGGGTGATCATTGTTACCTCTTTCCAGGTCTGGCATATAGGTATGTATGACTTGCATG rezultatas =>  Ešerichija marmotae - gramneigiamų, sporas nesudarančių, lazdelės formos bakterija. 
Detaliau galima rasti /3 Laboratorinis darbas/Peakai/Peakas 1 dalis 2.txt



Antras peak'as:

Sekos TTTGGTTCTTCCCGATAATCCTGGGATACACCGCGGGGAAACGCTTCGGCGGTAATCCATTTACTGCCATGGTGATTGGTGGAGCGTTAGTGCATCCATTAATTCTGACTGCTTTCGAGAACGGGCAAAAAGCGGATGCGCTGGGGCTGG rezultatas => Šigella sonnei - gramneigiama, lazdelės pavidalo, nemotorinė, sporą nesukelianti bakterija. 

Detaliau galima rasti /3 Laboratorinis darbas/Peakai/Peakas 2 dalis 1.txt

Sekos GGTGGTTTCGCCGGGCATGGTTGATGCGCACACCCATATTTCTGAACCGGGTCGTAGCCACTGGGAAGGTTATGAAACCGGTACTCGCGCAGCGGCAAAAGGTGGTATCACCACCATGATCGAAATGCCGCTCAACCAGCTGCCTGCAACG rezultatas => Ešerichija coli (žarninės lazdelės bakterija). 
Jau buvo Ešerichijos šeimos bakterija šitame mėginyje, galima daryt prielaidą, kad jų yra ir daugiau. 

Detaliau galima rasti /3 Laboratorinis darbas/Peakai/Peakas 2 dalis 2.txt



Trečias peak'as:

Sekos GACGTCCTGGACTGCGGCAACGCCGGAACCCTCATGCGCCTCCTCCTCGGCCTCCTCGCGGGCCAGGAGGGGCTTTTCGCCGTCCTCACCGGGGACGCCTCCTTGAGGCGCCGCCCCATGGGCCGGGTGGTGGCCCCCTTGAGGGCCATGG rezultatas => Thermus thermophilus - gramneigiama bakterija, naudojama daugelyje biotechnologinių taikymo sričių, įskaitant kaip pavyzdinį organizmą genetinei manipuliacijai, struktūrinę genomiką ir sistemų biologiją, bakterija yra ypač atspari karščiui. 

Detaliau galima rasti /3 Laboratorinis darbas/Peakai/Peakas 3 dalis 1.txt

Sekos GTCCACGGCGAAGGCCAGGTCCACCAGGTTGATCAGGAGGACCACCCGCCAGAACGTCTTCGCCGAGGCCTCGGGGAGGGGCTTCCCCTCGGGGTGGTTGCGGAAGTGCCGGACCATGAGGTAGACGAGGTAGAGCCCCCCGAGGACCTGG rezultatas => Thermus thermophilus - gramneigiama bakterija, naudojama daugelyje biotechnologinių taikymo sričių, įskaitant kaip pavyzdinį organizmą genetinei manipuliacijai, struktūrinę genomiką ir sistemų biologiją, bakterija yra ypač atspari karščiui. 
Tai nėra nuostabu, kad pasikartoja bakterija, nes šitame "peake" yra labai didelis kiekis C/G nukleotidų. 

Detaliau galima rasti /3 Laboratorinis darbas/Peakai/Peakas 3 dalis 2.txt

--------------------------------------------------------------------------------------------------------------------

5) Kokių rūšių bakterijų buvo mėginyje?

 - Bakterijų su dideliu atsparumu karštai temperatūrai (Thermo Thermophilus), tai įtakojo didelis C/G nukleotidų pasiskirstymas. 
 - Šiltakraujuose organizmuose randamų bakterijų (Pvz.: Šigella sonnei, Ešerichija coli). 

Buvo patikrinta tik dalis DNR sekų atitinkančių peak'us.

--------------------------------------------------------------------------------------------------------------------

# Antras laboratorinis darbas

Naudotos komandos užduoties 3 punkte:
perl script/sortgenome.pl --genomes-file my_data/NC_045512.fasta --sortedgenomes-file my_data/NC_045512.sort.fasta
perl script/sortgenome.pl --genomes-file my_data/MN514967.1.fasta --sortedgenomes-file my_data/MN514967.1.sort.fasta
perl script/sortgenome.pl --genomes-file my_data/combined_sequence.fasta --sortedgenomes-file my_data/combined_sequence.sort.fasta
perl script/sortgenome.pl --genomes-file my_data/complete_sequence.fasta --sortedgenomes-file my_data/complete_sequence.sort.fasta

./gclust -minlen 20 -both -nuc -threads 8 -ext 1 -sparse 2 -memiden 97 my_data/NC_045512.sort.fasta > my_data/NC_045512.sort.fasta.clustering.out
./gclust -minlen 20 -both -nuc -threads 8 -ext 1 -sparse 2 -memiden 97 my_data/MN514967.1.sort.fasta > my_data/MN514967.1.sort.fasta.clustering.out
./gclust -minlen 20 -both -nuc -threads 8 -ext 1 -sparse 2 -memiden 97 my_data/combined_sequence.sort.fasta > my_data/combined_sequence.sort.fasta.clustering.out
./gclust -minlen 20 -both -nuc -threads 8 -ext 1 -sparse 2 -memiden 97 my_data/complete_sequence.sort.fasta > my_data/complete_sequence.sort.fasta.clustering.out

seqkit grep -r -f my_data/NC_045512.fasta my_data/NC_045512.sort.fasta.clustering.out -o my_data/NC_045512.seqkit.extract.out
seqkit grep -r -f my_data/MN514967.1.fasta my_data/MN514967.1.sort.fasta.clustering.out -o my_data/MN514967.1.seqkit.extract.out
seqkit grep -r -f my_data/combined_sequence.fasta my_data/combined_sequence.sort.fasta.clustering.out -o my_data/combined_sequence.seqkit.extract.out
seqkit grep -r -f my_data/complete_sequence.fasta my_data/complete_sequence.sort.fasta.clustering.out -o my_data/complete_sequence.seqkit.extract.out

#Papildomai buvo pabandyta sugeneruoti genomus
make -f script/makefile_createreps gf=my_data/NC_045512.fasta clu=my_data/NC_045512.sort.fasta.clustering.out pgf=my_data/NC_045512.repgenomes.out CRG
make -f script/makefile_createreps gf=my_data/MN514967.1.fasta clu=my_data/MN514967.1.sort.fasta.clustering.out pgf=my_data/MN514967.1.repgenomes.out CRG
make -f script/makefile_createreps gf=my_data/combined_sequence.fasta clu=my_data/combined_sequence.sort.fasta.clustering.out pgf=my_data/combined_sequence.repgenomes.out CRG
make -f script/makefile_createreps gf=my_data/complete_sequence.fasta clu=my_data/complete_sequence.sort.fasta.clustering.out pgf=my_data/complete_sequence.repgenomes.out CRG

Naudotos komandos užduoties 4 punkte:
#https://bioinf.shenwei.me/seqkit/usage/#translate
seqkit translate --trim my_data/D3W8N4.against.viral.genomes.fasta > my_data/D3W8N4.against.viral.genomes.translated.fasta

#https://bioinf.shenwei.me/seqkit/usage/#seq
seqkit seq -m 800 my_data/D3W8N4.against.viral.genomes.translated.fasta > my_data/D3W8N4.against.viral.genomes.translated.filtered.fasta

#https://mafft.cbrc.jp/alignment/software/manual/manual.html
mafft --maxiterate 1000 --localpair my_data/D3W8N4.against.viral.genomes.translated.filtered.fasta > my_data/D3W8N4.aligned.fasta 

#https://www.biostars.org/p/170941/
#https://www.theunixschool.com/2014/08/sed-examples-remove-delete-chars-from-line-file.html
perl -i.bak -pe 's/\h+$//' my_data/D3W8N4.aligned.fasta
perl -i.bak -pe 's/://' my_data/D3W8N4.aligned.fasta

#http://www.microbesonline.org/fasttree/
fasttree -gamma my_data/D3W8N4.aligned.fasta > my_data/tree.txt

#http://etetoolkit.org/treeview/
Paste my_data/tree.txt

Naudotos komandos užduoties 5 punkte:

#http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html    (set_outgroup)

Google colab link:
https://colab.research.google.com/drive/16TCoEgAs2bJuLbgUVQAfLmfggV_0JF8o#scrollTo=pCWclyFTXjtE
--------------------------------------------------------------------------------------------------------------------

Interpretation.....how did the Covid-19 evolve, what path through hosts was taken? Would it be
different interpretation if out-group is not used? What about Urbani SARS origin? Is the Palm Civet
origin evident?

--------------------------------------------------------------------------------------------------------------------

1) Interpretation.....how did the Covid-19 evolve, what path through hosts was taken?

Kaip matome iš paskutinio gauto medžio, Covid-19 evoliuconavo ne iš kupranugarių covid viruso. O perėjo per šiuos virusus KF293666.1 / KF293664.1 / MW039392.1 / MT085168.1 / MT085175.1 / MT085172.1 / MZ064540.1 / MZ064531.1 / LC469308.1 / NC_034440.1. Bet iš kitų šaltinių radau, kad kupranugarių virusas yra giminingas covid-19. Tikėtina, kad skaičiavimai buvo atlikti netaisyklingai, įsivėlė klaida.

--------------------------------------------------------------------------------------------------------------------

2) Would it be different interpretation if out-group is not used?

Taip, kitu atveju matysis (būtų traktuojami) skirtingai giministės ryšiai in/out-group'sų.

--------------------------------------------------------------------------------------------------------------------
3) What about Urbani SARS origin?

Iš dabartinio medžio nėra matoma Urbani SARS kilmė. Tačiau iš kitų šaltinių radau, kad Urbani SARS kilo nuo šikšnosparnių ir buvo perneštas žmogui per laukines kates.

--------------------------------------------------------------------------------------------------------------------

4) Is the Palm Civet origin evident?

Iš dabartinio medžio nėra aišku. Negalima atsakyti.

--------------------------------------------------------------------------------------------------------------------

Atsakymai gali būti neteisingi arba dalinai teisingi dėl atsiradusių klaidų skaičiavimuose ar genomų duombazei grąžinus klaidingus duomenis. Visi tarpiniai sprendimai, duomenys bei rezultatai yra sukelti prie atitinkamų punktų.

# Pirmas laboratorinis darbas

Visi dokumentai ir kodas susijęs su šia užduotimi yra "1 Laboratorinis darbas" aplanke. 

1 Laboratorinis darbas/inputs/ aplanke laikomi fasta dokumentai ir užduoties aprašymas.

1 Laboratorinis darbas/outputs/ aplanke laikomi kodonų ir dikodonų medžiai, bei programos išvesties tekstas.

1 Laboratorinis darbas/libs/ aplanke laikomi dependencies jar'ai.

1 Laboratorinis darbas/src/ aplanke laikomas programos veikimo kodas.

Metodas "List<DNASequence> findCodonSequences(List<DNASequence> seq)" yra naudojamas atrasti visus reikalingus sekų
trinarius, kurie prasideda pradžios kodonu "ATG" ir baigiasi vienu iš pabaigos kodonu "TAA, "TAG" ar "TGA".

Metodas "void printSequenceList(String text, List<DNASequence> seq)" yra naudojamas atspausdinti kodonų seką.

Metodas "DNASequence findLongestCodonSequence(List<DNASequence> seq)" yra naudojamas atrasti ilgiausią kodonų seką (2
užduotis).

Metodas "void findCodonSequence(List<DNASequence> seq)" yra naudojamas įvertinti kodonų dažnius (4 užduotis).

Metodas "void findDiCodonSequence(List<DNASequence> seq)" yra naudojamas įvertinti dikodonų dažnius (4 užduotis).

Metodas "List<DNASequence> removeUselessCodons(List<DNASequence> seq)" yra naudojamas pašalinti visas nereikalingas
kodonų sekas, kurios trumesnės nei 100, kad turėtume realių DNA sekų rezultatus.

Metodas "void calculateCodonMatrix()" yra naudojamas apskaičiuoti ir atspausdinti kodonų matricą<br/>
Buvo naudojama ši formulė apskaičiuoti kodonų matricą<br/>

F1 - pirmas failas, F2 - antras failas<br/>
(ATG[F1] - ATG[F2])^2 = x1<br/>
(TAA[F1] - TAA[F2])^2 = x2<br/>
(TAG[F1] - TAG[F2])^2 = x3<br/>
(TGA[F1] - TGA[F2])^2 = x4<br/>

x1 + x2 + x3 + x4 - suma, kuria įstatom į matricą.<br/>

Philip formatu (matrica) atrodytų:

8<br/>
b1 X, skaičius(b2), skaičius(b3), skaičius(b4), skaičius(m1), skaičius(m2), skaičius(m3), skaičius (m4)<br/>
b2 skaičius(b1), X, skaičius(b3), skaičius(b4), skaičius(m1), skaičius(m2), skaičius(m3), skaičius (m4)<br/>
b3 ...<br/>
b4 ...<br/>
m1 ...<br/>
m2 ...<br/>
m3 ...<br/>
m4 ...<br/>

Metodas "void calculateDiCodonMatrix()" atitinkamai skaičiuoja dikodonų dažnį. Naudojama tokia pati formulė kaip ir
skaičiuojant kodonus. Dikodonuos skaičiuojami juos transliuojant iš DNA sekos. Atrasti trinariai pavaizduojami vienu
simboliu panaudojant TranscriptionEngine.translate() (Biojava) metodą. Tada transliuoti simboliai yra sujungiami po du
ir randamas dažnis tokių simbolio porų. Toliau atliekamas skaičiavimas (toks pat kaip kodonams skaičiuoti) ir sudedama į
simbolių matricą:<br>
(XX[F1] - XX[F2])^2

Metodas "void execute(String fileName, int number)" yra naudojamas įvykdyti 1 - 5 užduotis.

Pagalbiniai metodai atlikti skaičiavimams:

Metodas "int getSpecificValueFromDictionary(int z1, int z3)" yra naudojamas grąžinti diCodonų sekos atitinkamam nariui.

Utils klasė:

Metodas "DNASequence readFastaFile(String fileName)" yra naudojamas nuskaityti Fasta tipo dokumento turinį.

Metodas "List<DNASequence> splitStringToFrames(int position, String data, int split)" yra naudojamas suskaidyti DNA
seką (String)
atitinkamais dydžiais, nuo tam tikros pozicijos.

Metodas "Map<String, Integer> countOccurrences(List<String> splitSeq)" yra naudojamas suskaičiuoti sekų pasikartojimus.

Metodas "<K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map)" yra naudojamas surūšiuoti kolekcijai
pagal reikšmes.

Metodas "List<String> dnaSequenceToString(List<DNASequence> seq)" yra naudojamas paversti sąrašo elementus iš
DNASequence į String formatą.

Metodas "TranscriptionEngine getTranscriptionEngine()" yra naudojamas gauti TranscriptionEngine objektą.

# Pirmojo laboratorinio darbo analizė

Programos veikimo rezultatai:

Gauta kodonų matricą:

B1 B2 B3 B4 M1 M2 M3 M4<br/>
0.000000 0.000068 0.000022 0.000077 0.000003 0.000011 0.000084 0.000008 <br/>
0.000068 0.000000 0.000014 0.000282 0.000067 0.000031 0.000011 0.000070 <br/>
0.000022 0.000014 0.000000 0.000179 0.000022 0.000006 0.000022 0.000031 <br/>
0.000077 0.000282 0.000179 0.000000 0.000078 0.000129 0.000320 0.000078 <br/>
0.000003 0.000067 0.000022 0.000078 0.000000 0.000009 0.000085 0.000009 <br/>
0.000011 0.000031 0.000006 0.000129 0.000009 0.000000 0.000050 0.000011 <br/>
0.000084 0.000011 0.000022 0.000320 0.000085 0.000050 0.000000 0.000104 <br/>
0.000008 0.000070 0.000031 0.000078 0.000009 0.000011 0.000104 0.000000 <br/>

Konvertuota kodonų matrica į phylip formatą:

8 <br/>
B1 0.000068 0.000022 0.000077 0.000003 0.000011 0.000084 0.000008 <br/>
B2 0.000068 0.000014 0.000282 0.000067 0.000031 0.000011 0.000070 <br/>
B3 0.000022 0.000014 0.000179 0.000022 0.000006 0.000022 0.000031 <br/>
B4 0.000077 0.000282 0.000179 0.000078 0.000129 0.000320 0.000078 <br/>
M1 0.000003 0.000067 0.000022 0.000078 0.000009 0.000085 0.000009 <br/>
M2 0.000011 0.000031 0.000006 0.000129 0.000009 0.000050 0.000011 <br/>
M3 0.000084 0.000011 0.000022 0.000320 0.000085 0.000050 0.000104 <br/>
M4 0.000008 0.000070 0.000031 0.000078 0.000009 0.000011 0.000104 <br/>

Gauta dikodonų matricą:

B1 B2 B3 B4 M1 M2 M3 M4<br/>
0.000000 0.001408 0.001426 0.001322 0.001802 0.001979 0.001539 0.002768 <br/>
0.001408 0.000000 0.000347 0.001660 0.000460 0.000387 0.000766 0.000908 <br/>
0.001426 0.000347 0.000000 0.001797 0.000408 0.000722 0.000820 0.001409 <br/>
0.001322 0.001660 0.001797 0.000000 0.002124 0.002151 0.002062 0.002822 <br/>
0.001802 0.000460 0.000408 0.002124 0.000000 0.000594 0.000887 0.001198 <br/>
0.001979 0.000387 0.000722 0.002151 0.000594 0.000000 0.001464 0.000344 <br/>
0.001539 0.000766 0.000820 0.002062 0.000887 0.001464 0.000000 0.002201 <br/>
0.002768 0.000908 0.001409 0.002822 0.001198 0.000344 0.002201 0.000000 <br/>

Konvertuota dikodonų matrica į phylip formatą:

8 <br/>
B1 0.001408 0.001426 0.001322 0.001802 0.001979 0.001539 0.002768 <br/>
B2 0.001408 0.000347 0.001660 0.000460 0.000387 0.000766 0.000908 <br/>
B3 0.001426 0.000347 0.001797 0.000408 0.000722 0.000820 0.001409 <br/>
B4 0.001322 0.001660 0.001797 0.002124 0.002151 0.002062 0.002822 <br/>
M1 0.001802 0.000460 0.000408 0.002124 0.000594 0.000887 0.001198 <br/>
M2 0.001979 0.000387 0.000722 0.002151 0.000594 0.001464 0.000344 <br/>
M3 0.001539 0.000766 0.000820 0.002062 0.000887 0.001464 0.002201 <br/>
M4 0.002768 0.000908 0.001409 0.002822 0.001198 0.000344 0.002201 <br/>

Kodonų(CodonTree.svg) ir dikodonų(DiCodonTree.svg) medžių nuotraukos yra sukeltos į "BioInformatika\outputs" aplanką

Dikodonų medyje (DiCodonTree.svg) matoma, klasterizacija į 4 dalis, dalys yra žymiai didesnio atstumo, negu kodonų
dažnių medyje. Kodonų medyje (CodonTree.svg) matoma, klasterizacija į 4 dalis, pačios atšakos yra minimalaus atstumo
tarpusavyje. Minėtų klasterizacijų 4 dalys yra B4 ir B3, B1 ir M4, M1 ir M3, B2 ir M2.

Tarp dviejų medžių matomas didelis skirtumas - dažnis tarp kodonų ir dikodonų skiriasi. Tai gali būti dėl skaičiavimų
netikslumų. Taip pat klasteriuose yra pastebima, kad kai kur bakterijų ir žinduolių virusų dažniai sutampa. Labiausiai,
varijuoja B4 virusas tiek kodonų, tiek dikodonų medyje.
