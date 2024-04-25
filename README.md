# Detailed Sequences for Trees Unblemished

> Detailed Sequences for Trees Unblemished (DSTU) simplifies phylogenetic tree creation in microbiology and virology by facilitating sequence downloads from NCBI GenBank using accession numbers. It also reinstates organism names in trees constructed with IQ-TREE, retrieves host information about microorganisms, and prepares annotation datasets for iTOL.

|The Good|The Bad|The Ugly|
|--------|-------|--------|
|<img src="https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished/blob/main/imgs/second%20tree.png" width="100%">|<img src="https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished/blob/main/imgs/third%20tree.png" width="100%">|<img src="https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished/blob/main/imgs/first%20tree.png" width="100%">|

DSTU allows easy and simple annotation of phylogenetic trees. The best tree contains information about the hosts from which the virus was isolated and the full names of the viruses. The so-so tree contains the same information, but is colour annotated with randomly generated colours. The worst tree contains only accession numbers on its leaves.

## Features
### Sequence Downloading
- Facilitates the retrieval of sequences from NCBI GenBank using specified accession numbers.
### Organism Name Reintegration
- Enhances IQ-TREE constructed phylogenetic trees by replacing accession numbers with the corresponding organism names for clarity and context.
### Host Information Retrieval
- Gathers comprehensive host data for each microorganism, including the host's phylogenetic order.
### Annotation Dataset Preparation for iTOL
- Utilizes the collected host information to prepare detailed annotation datasets, optimizing visualization in iTOL.

## Installation

```bash
git clone https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished.git && cd Detailed-Sequences-for-Trees-Unblemished
```

```bash
pip install biopython
```

## Usage Guide
- Demonstrational python notebook is available in `demo.ipynb` file
- Demonstrational data is available in `demo_data` folder

Demonstrational data is based on the recent paper about identifying novel hantavirus in bats

### Import modules

**_Input_**

```python
from DSTU import *
```

Further we will need this short function:

```python
def print_first_5_lines(filename):
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if i >= 5:
                break
            print(line.strip())
```

This function has no sense regarding DSTU tool. It will be used to demonstrate file contents.

To use DSTU user must have `accession numbers.txt` file. It must look like this:

**_Input_**

```python
print_first_5_lines('demo_data/accession_numbers.txt')
```

**_Output_**

```
NC_034519
NC_055636
NC_005225
NC_038939
NC_038529
```

### Example of using the `get_sequences` function

**_Input_**

```python
get_sequences('iljapopov17@gmail.com', 'demo_data/accession_numbers.txt', 'genbank_sequences')
```

Input:
1. email
2. input txt file with the list of accession numbers
3. output file

**_Output_**

```
Downloaded: NC_034519
Downloaded: NC_055636
Downloaded: NC_005225
Downloaded: NC_038939
Downloaded: NC_038529
Downloaded: NC_034467
Downloaded: NC_034553
Downloaded: NC_003468
Downloaded: NC_034515
Downloaded: NC_038299
Downloaded: NC_077671
Downloaded: NC_034403
Downloaded: NC_038695
Downloaded: NC_034556
Downloaded: LC553715
Downloaded: NC_005238
Downloaded: NC_005235
Downloaded: NC_034517
Downloaded: NC_006435
Downloaded: NC_005222
Downloaded: NC_055147
Downloaded: NC_034560
Downloaded: NC_034485
Downloaded: NC_034407
Downloaded: NC_034399
Downloaded: NC_034402
Downloaded: MG663536
Downloaded: NC_038515
Downloaded: NC_078262
Downloaded: NC_055632
Downloaded: NC_034401
Downloaded: OR684449
Downloaded: FJ593498
Downloaded: KX512433
Downloaded: NC_078485
Downloaded: KX779126
Downloaded: NC_034564
Downloaded: NC_010707
Downloaded: NC_055170
All downloads completed.
```

### What's next? Tree construction.

Below used tree is based on demonstrational data. Tree construction is made in `tree_construction.ipynb` notebook.<br>
The most important file is `demo_data/tree_ufb.treefile`. Upload it to iTOL for visualization.

<img src="https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished/blob/main/imgs/Reference%20tree.jpg" width="75%">

_Figure 1. Reference tree from the original paper_

<img src="https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished/blob/main/imgs/first%20tree.png" width="75%">

_Figure 2. Naked phylogenetic tree_

This tree is naked.<br>
There is no:<br>
1. Annotation of the organisms name. There are only accession numbers that cannot say anything.
2. The tree demonstrates phylogenetic relationships between different viruses. But there is no information about host organisms of that viruses.

It is worth mentioning that the trees are literally identical. (Bootstrap values are even better in my variant).

### Example of returning organisms names to the tree with `get_organisms` function

**_Input_**

```python
get_organisms('iljapopov17@gmail.com', 'demo_data/accession_numbers.txt', 'demo_data/accession_organism.txt')
```

Input:
1. email
2. input txt file with the list of accession numbers
3. output file

**_Output_**

```
The request has been fulfilled.
File saved to demo_data/accession_organism.txt
```

**_Input_**

```python
print_first_5_lines('demo_data/accession_organism.txt')
```

**_Output_**

```
NC_034519.1 Orthohantavirus khabarovskense
NC_055636.1 Orthohantavirus tatenalense
NC_005225.1 Orthohantavirus puumalaense
NC_038939.1 Orthohantavirus prospectense
NC_038529.1 Eothenomys miletus hantavirus LX309
```

**_Input_**

```python
update_tree('demo_data/accession_organism.txt', 'demo_data/tree_ufb.treefile', 'demo_data/modified_tree.treefile')
```

Input:
1. input txt file with the list of accession numbers and organisms names
2. input tree file
3. output modified tree

**_Output_**

```
The request has been fulfilled.
File saved to demo_data/accession_organism.txt
```

**_Input_**

```python
print_first_5_lines('demo_data/tree_ufb.treefile')
```

I know it was a stupid idea to read first 5 lines in a tree file... Anyway. It is made in demonstrational purposes.

**_Output_**

```
(FJ593498.1:0.1240225441,KX512433.1:0.1580233515,((((((KX779126.1:0.1801341369,NC_034564.1:0.1518834757)100:0.2690724126,NC_010707.1:0.4026159852)100:0.5357342048,NC_055170.1:3.2821188681)96:0.1993926731,(((((((LC553715.1:0.2424396410,NC_034556.1:0.2425091493)100:0.1540926638,NC_005238.1:0.2987153355)100:0.0815585291,((NC_005222.1:0.1745023294,NC_006435.1:0.1329576555)100:0.2619343862,(NC_005235.1:0.3091225291,NC_034517.1:0.3426538757)100:0.0890215926)59:0.0409665566)100:0.2464858691,NC_055147.1:0.5009202874)69:0.0579848758,((NC_034399.1:0.4574934455,NC_034407.1:0.4201827666)100:0.2133110745,(NC_034485.1:0.3554924125,NC_034560.1:0.3958031671)100:0.1134597575)100:0.0983251229)88:0.0539862361,NC_034402.1:0.6194957047)100:0.2179091508,(((((NC_003468.2:0.3220309131,NC_034553.1:0.3217768427)100:0.0967750566,(NC_034515.1:0.3420020277,NC_038299.1:0.3578938480)78:0.0604905717)100:0.0681808060,(NC_034403.1:0.4057149461,NC_077671.1:0.3295415521)96:0.0808310506)100:0.1515347499,NC_038695.1:0.6146030315)75:0.0529123693,(((NC_005225.1:0.3178994625,(NC_034519.1:0.2903408237,NC_055636.1:0.2951103060)96:0.0707049689)100:0.1162566928,NC_038939.1:0.4860686808)100:0.0974816090,(NC_034467.1:0.3408088379,NC_038529.1:0.3214413064)100:0.1876647016)100:0.0906674433)100:0.3112288106)100:0.2995543026)97:0.1136359007,NC_078485.1:1.2137610889)49:0.0697049196,(((MG663536.1:0.4927348232,NC_038515.1:0.3837609395)94:0.0895431598,NC_078262.1:0.4767046102)100:0.2182159381,((NC_034401.1:0.5482148765,NC_055632.1:0.5333969980)100:0.2727779310,OR684449.1:0.6549294470)90:0.1135643862)55:0.0661132415)100:0.9075896851);
```

Usual treefile contains only accession numbers. They cannot say anything.

**_Input_**

```python
print_first_5_lines('demo_data/modified_tree.treefile')
```

**_Output_**

```
(FJ593498.1 Nova virus:0.1240225441,KX512433.1 Nova virus:0.1580233515,((((((KX779126.1 Imjin virus:0.1801341369,NC_034564.1 Imjin virus:0.1518834757)100:0.2690724126,NC_010707.1 Thottapalayam virus:0.4026159852)100:0.5357342048,NC_055170.1 Hainan oriental leaf-toed gecko hantavirus:3.2821188681)96:0.1993926731,(((((((LC553715.1 Orthohantavirus thailandense:0.2424396410,NC_034556.1 Anjozorobe virus:0.2425091493)100:0.1540926638,NC_005238.1 Orthohantavirus seoulense:0.2987153355)100:0.0815585291,((NC_005222.1 Orthohantavirus hantanense:0.1745023294,NC_006435.1 Hantavirus Z10:0.1329576555)100:0.2619343862,(NC_005235.1 Orthohantavirus dobravaense:0.3091225291,NC_034517.1 Orthohantavirus sangassouense:0.3426538757)100:0.0890215926)59:0.0409665566)100:0.2464858691,NC_055147.1 Tigray virus:0.5009202874)69:0.0579848758,((NC_034399.1 Jeju virus:0.4574934455,NC_034407.1 Bowe virus:0.4201827666)100:0.2133110745,(NC_034485.1 Orthohantavirus caobangense:0.3554924125,NC_034560.1 Kenkeme virus:0.3958031671)100:0.1134597575)100:0.0983251229)88:0.0539862361,NC_034402.1 Bruges virus:0.6194957047)100:0.2179091508,(((((NC_003468.2 Orthohantavirus andesense:0.3220309131,NC_034553.1 Maporal virus:0.3217768427)100:0.0967750566,(NC_034515.1 Orthohantavirus delgaditoense:0.3420020277,NC_038299.1 Orthohantavirus bayoui:0.3578938480)78:0.0604905717)100:0.0681808060,(NC_034403.1 Orthohantavirus montanoense:0.4057149461,NC_077671.1 Orthohantavirus sinnombreense:0.3295415521)96:0.0808310506)100:0.1515347499,NC_038695.1 Rockport virus:0.6146030315)75:0.0529123693,(((NC_005225.1 Orthohantavirus puumalaense:0.3178994625,(NC_034519.1 Orthohantavirus khabarovskense:0.2903408237,NC_055636.1 Orthohantavirus tatenalense:0.2951103060)96:0.0707049689)100:0.1162566928,NC_038939.1 Orthohantavirus prospectense:0.4860686808)100:0.0974816090,(NC_034467.1 Fugong virus:0.3408088379,NC_038529.1 Eothenomys miletus hantavirus LX309:0.3214413064)100:0.1876647016)100:0.0906674433)100:0.3112288106)100:0.2995543026)97:0.1136359007,NC_078485.1 Lena virus:1.2137610889)49:0.0697049196,(((MG663536.1 Dakrong virus:0.4927348232,NC_038515.1 Laibin virus:0.3837609395)94:0.0895431598,NC_078262.1 Xuan son virus:0.4767046102)100:0.2182159381,((NC_034401.1 Quezon virus:0.5482148765,NC_055632.1 Orthohantavirus robinaense:0.5333969980)100:0.2727779310,OR684449.1 Buritiense virus:0.6549294470)90:0.1135643862)55:0.0661132415)100:0.9075896851);
```

Modified treefile contains accession numbers and organisms names. It makes more sense.

### Example of fetching hosts info with `get_hosts` function

**_Input_**

```python
get_hosts('iljapopov17@gmail.com', 'demo_data/accession_numbers.txt', 'demo_data/accession_host.txt')
```

Input:
1. email
2. input txt file with the list of accession numbers
3. output file

**_Output_**

```
The request has been fulfilled.
File saved to demo_data/accession_host.txt
```

**_Input_**

```python
print_first_5_lines('demo_data/accession_host.txt')
```

**_Output_**

```
NC_034519.1 Microtus maximowiczii
NC_055636.1 Microtus agrestis
NC_005225.1 ND
NC_038939.1 Microtus pennsylvanicus
NC_038529.1 Eothenomys miletus
```

**_Input_**

```python
get_hosts_orders('iljapopov17@gmail.com', 'demo_data/accession_host.txt', 'demo_data/accession_order.txt')
```

Input:
1. email
2. input txt file with the list of accession numbers and hosts
3. output file

**_Output_**

```
The request has been fulfilled.
Please do not forget to edit the file manually.
The query to NCBI database from this function is pretty difficult.
Sometimes this function prints:
"Error - HTTP Error 400: Bad Request" in case of bad connection or
"Note - False record" in case there is no record about the host organism.
```

**_Input_**

```python
print_first_5_lines('demo_data/accession_order.txt')
```

**_Output_**

```
NC_034519.1	Rodentia
NC_055636.1	Rodentia
NC_005225.1	ND
NC_038939.1	Rodentia
NC_038529.1	Rodentia
```

### Example of preparing info for iTOL

**_Input_**

```python
unique_orders = get_unique_orders("demo_data/accession_order.txt")
print(unique_orders)
```

**_Output_**

```
['Rodentia', 'ND', 'Eulipotyphla', 'Chiroptera', 'Squamata']
```

**_Input_**

```python
color_map = set_color_map("demo_data/accession_order.txt")
print(color_map)
```

Interactive window will open and will ask to set HEX codes for each unique order

**_Output_**

```
{'Rodentia': '#0ca20c', 'ND': '#ffffff', 'Eulipotyphla': '#0078ff', 'Chiroptera': '#000000', 'Squamata': '#ffa500'}
```

### Example of creating annotation dataset for iTOL

#### Using the manually adjusted color map

**_Input_**

```python
# Example usage:
# To use with a user-defined color map:
get_itol_dataset("demo_data/accession_organism.txt", "demo_data/accession_order.txt", "demo_data/dataset_for_iTOL.txt", color_map)
```

Input:
1. input txt file with the list of accession numbers and organisms names
2. input txt file with the list of accession numbers and phylogenetic order of microorganism host
3. output file
4. manually created color map

**_Output_**

```
Colors were set by the user.
The request has been fulfilled.
```

**_Input_**

```python
print_first_5_lines('demo_data/dataset_for_iTOL.txt')
```

**_Output_**

```
DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	Host Group Colors
DATA
NC_034519.1 Orthohantavirus khabarovskense	#0ca20c	Rodentia
```

#### Next steps

1. Visit iTOL
2. Upload `demo_data/modified_tree.treefile` file as the tree
3. Upload `demo_data/dataset_for_iTOL.txt` as the annotation dataset

<img src="https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished/blob/main/imgs/second%20tree.png" width="75%">

https://itol.embl.de/export/95255816295731714078181

_Fig 2. Second tree. With annotation info containing organisms names and manually adjusted colors indicating hosts phylogenetic order_

This is the best tree easily made with DSTU software

|Original|DSTU|
|--------|-------|
|<img src="https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished/blob/main/imgs/Reference%20tree.jpg">|<img src="https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished/blob/main/imgs/second%20tree.png">|

It can be seen that in original version authors did annotation manually and they made some mistakes in hosts annotation. DSTU software did not make this mistakes.

#### Using randomly generated color map

**_Input_**

```python
get_itol_dataset("demo_data/accession_organism.txt", "demo_data/accession_order.txt", "demo_data/dataset_for_iTOL_2.txt")
```

Input:
1. input txt file with the list of accession numbers and organisms names
2. input txt file with the list of accession numbers and phylogenetic order of microorganism host
3. output file

**_Output_**

```
Colors were not set, they were generated randomly.
The request has been fulfilled.
```

**_Input_**

```python
print_first_5_lines('demo_data/dataset_for_iTOL_2.txt')
```

**_Output_**

```
DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	Host Group Colors
DATA
NC_034519.1 Orthohantavirus khabarovskense	#e31342	Rodentia
```

#### Next steps

1. Visit iTOL
2. Upload `demo_data/modified_tree.treefile` file as the tree
3. Upload `demo_data/dataset_for_iTOL_2.txt` as the annotation dataset

<img src="https://github.com/iliapopov17/Detailed-Sequences-for-Trees-Unblemished/blob/main/imgs/third%20tree.png" width="75%">

https://itol.embl.de/export/95255816343181714079366

_Fig 3. Third tree. With annotation info containing organisms names and randomly generated colors indicating hosts phylogenetic order_

In this case random generation played a bad joke! Almost every color is the same. It will be much more convenient to adjust color map manually.
