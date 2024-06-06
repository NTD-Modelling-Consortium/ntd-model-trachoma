import unittest

import numpy as np

import trachoma.trachoma_functions as tf

class TestTrachomaFunctions(unittest.TestCase):

    def test_initialize_dataframe(self):
        # Arrange
        outputYear = range(2016, 2017) # one year, 2016
        results = [(1, 2, 3, 4), (5, 6, 7, 8)] # two simulations with dummy data

        # Act
        df = tf.initialize_dataframe(outputYear, results)

        # Assert on the structure of the dataframe
        column_names = ["Time", "age_start", "age_end", "measure", 4, 5]
        self.assertListEqual(list(df.columns), column_names)
        # Assert on columns of the dataframe
        for column in column_names:
            self.assertListEqual(list(np.zeros(240)), list(df.loc[:, column]))

    def test_update_measure_for_index(self):
        # Arrange
        year = 2016
        outputYear = range(year, year + 1)  # one year, 2016
        results = [(1, 2, 3, 4), (5, 6, 7, 8)]  # two simulations with dummy data
        df = tf.initialize_dataframe(outputYear, results)
        max_age = 240 # why is this called max_age?
        measure = tf.Measure.HEAVY_INFECTIONS
        infs = np.ones(max_age) # array of dummy infection data
        nums = np.ones(max_age) # what is nums? dummy for now

        # Act, for index = 0, i = 0
        index = 0
        i = 0
        tf.update_measure_for_index(df, measure, index, max_age, infs, nums, i, year)

        # Assert that column 4 is populated with infection data
        self.assertListEqual(list(infs), list(df.loc[:, 4]))

        # Act, for index = 0, i = 1
        index = 0
        i = 1
        tf.update_measure_for_index(df, measure, index, max_age, infs, nums, i, year)

        # Assert that column 4 and 5 are populated with infection data
        self.assertListEqual(list(infs), list(df.loc[:, 4]))
        self.assertListEqual(list(infs), list(df.loc[:, 5]))
