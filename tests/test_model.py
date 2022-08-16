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

    for i, vals in enumerate(sim_outputs):
        expected_final_pop_state = numpy.loadtxt(f'tests/reference_output/final_pop_state_{i}.txt', delimiter=',').reshape(1000, 3)
        expected_yearly_infect_prev = numpy.loadtxt(f'tests/reference_output/yearly_infect_prev_{i}.txt', delimiter=',').reshape(2704, 2)
        expected_yearly_threshold_infs = numpy.loadtxt(f'tests/reference_output/yearly_threshold_infs_{i}.txt', delimiter=',').reshape(2705, 60)

        expected = {}
        expected['IndI'] = expected_final_pop_state[:,0]
        expected['IndD'] = expected_final_pop_state[:,1]
        expected['No_Inf'] = expected_final_pop_state[:,2]
        expected['True_Prev_Disease_children_1_9'] = expected_yearly_infect_prev[:,0]
        expected['True_Infections_Disease_children_1_9'] = expected_yearly_infect_prev[:,1]
        expected['Yearly_threshold_infs'] = expected_yearly_threshold_infs

        numpy.testing.assert_array_equal(vals['IndI'], expected['IndI'])
        numpy.testing.assert_array_equal(vals['IndD'], expected['IndD'])
        numpy.testing.assert_array_equal(vals['No_Inf'], expected['No_Inf'])
        numpy.testing.assert_allclose(vals['True_Prev_Disease_children_1_9'], expected['True_Prev_Disease_children_1_9'], atol=5.e-7, rtol=2.e-4)
        numpy.testing.assert_allclose(vals['True_Infections_Disease_children_1_9'], expected['True_Infections_Disease_children_1_9'], atol=5.e-7, rtol=2.e-4)
        numpy.testing.assert_allclose(vals['Yearly_threshold_infs'], expected['Yearly_threshold_infs'], atol=5.e-7, rtol=2.e-4)
