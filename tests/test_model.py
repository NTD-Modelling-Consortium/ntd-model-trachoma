from trachoma import Trachoma_Simulation

def main():
    Trachoma_Simulation("tests/data/beta_values.csv",
                        "tests/data/mda_input_2008_2017.csv",
                        "prevalence_job.csv",
                        "infection_job.csv",
                        SaveOutput = False,
                        OutSimFilePath = None,
                        InSimFilePath = None,
                        logger = None)

if __name__ == "__main__":
    main()
