import pytest
from InfluenceCalculator import InfluenceCalculator
import pandas as pd
import os
from tqdm import tqdm


def test_run_multiple_seeds():

    print("\nInitializing InfluenceCalculator object")
    ic = InfluenceCalculator('tests/toy_network_example.sqlite')
    print("InfluenceCalculator object created")

    # Define seed
    meta_column = 'seed_01'
    seed_categories = pd.unique(ic.meta[meta_column])
    seed_categories = [element for element in seed_categories
                       if element not in [None, '']]

    # Get neuron ids to inhibit (sensory neurons in this case)
    silenced_neurons = ic.meta[
        ic.meta['super_class'].isin(['sensory',
                                     'ascending_sensory'])].root_id
    
    # Create the influence folder if it doesn't exist
    folder_name = "Influence"
    os.makedirs(folder_name, exist_ok=True)

    for seed_category in tqdm(seed_categories,
                              desc="Calculating influence scores"):
        # Get seed neuron ids
        seed_ids = ic.meta[ic.meta[meta_column] == seed_category].root_id 

        influence_df = None
        influence_df = ic.calculate_influence(seed_ids, silenced_neurons)

        # Save the DataFrame to a CSV file
        influence_df.to_csv(f"{folder_name}/{seed_category}_influence.csv",
                            index=False)
        
    print("Influence scores saved in: ",
          os.path.join(os.getcwd(), folder_name))


if __name__ == "__main__":
    """
    Example use:
    1. Create InfluenceCalculator object using file path where dataset is saved
    4. Define seed categories
    5. Define silenced neuron ids
    6. Define seed neuron ids (per category)
    7. Build influence dataframe with influence scores (per category)
    8. Save dataframe (per category)
    9. Loop over all seed catgories
    """
    pytest.main()
    
