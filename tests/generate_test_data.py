from trachoma import Trachoma_Simulation
import numpy.testing

def main():
    Trachoma_Simulation("tests/data/beta_values.csv",
                        "tests/data/mda_input_2008_2017.csv",
                        "prevalence_job.csv",
                        "infection_job.csv",
                        SaveOutput = True,
                        OutSimFilePath = "results/out.p",
                        InSimFilePath = None,
                        logger = None)

if __name__ == "__main__":
    import pickle
    main()
    with open("results/out.p", "rb") as f:
        sim_outputs = pickle.load(f)
    final_pop_state = numpy.empty((1000, 3))
    yearly_infect_prev = numpy.empty((2704, 2))
    for i, vals in enumerate(sim_outputs):
        final_pop_state[:,0] = vals['IndI']
        final_pop_state[:,1] = vals['IndD']
        final_pop_state[:,2] = vals['No_Inf']
        yearly_infect_prev[:,0] = vals['True_Prev_Disease_children_1_9']
        yearly_infect_prev[:,1] = vals['True_Infections_Disease_children_1_9']
        with open(f"results/final_pop_state_{i}.txt", "w") as f:
            final_pop_state.tofile(f, sep=',', format='%i')
        with open(f"results/yearly_infect_prev_{i}.txt", "w") as f:
            yearly_infect_prev.tofile(f, sep=',', format='%f')
