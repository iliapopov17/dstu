from Bio import Entrez


def read_accession_file(filename):
    """Read accession numbers from a file."""
    with open(filename, "r") as file:
        accession_numbers = [line.strip() for line in file.readlines()]
    return accession_numbers


def fetch_organism_names(email, accession_numbers):
    """Fetch organism names and versions from NCBI's Entrez for given accession numbers."""
    Entrez.email = email
    organism_dict = {}
    for accession in accession_numbers:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["GBSeq_organism"]
            accession_version = records[0][
                "GBSeq_accession-version"
            ]  # Fetch the accession version
            organism_dict[accession_version] = (
                organism_name  # Use accession version as key
            )
            handle.close()
        except IndexError:
            print(f"No data found for {accession}")
        except KeyError:
            print(f"Organism name or accession version missing for {accession}")
        except Exception as e:
            print(f"Error retrieving data for {accession}: {e}")
            handle.close()
    return organism_dict


def get_organisms(email, input_filename, output_filename):
    """Main function to get organisms with version and write them to a file."""
    accession_numbers = read_accession_file(input_filename)
    organism_dict = fetch_organism_names(email, accession_numbers)

    with open(output_filename, "w") as file:
        for accession_version, organism in organism_dict.items():
            file.write(f"{accession_version} {organism}\n")
