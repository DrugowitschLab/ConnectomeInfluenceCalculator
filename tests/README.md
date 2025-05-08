Running the tests additionally requires the `pytest` and `tqdm` packages installed.
To run the tests in this folder, navigate to the main `ConnectomeInfluenceCalculator` folder, and run
```sh
pytest -s tests
```
in the terminal.
This will load a simplified example connectome from `toy_network_example.sqlite`, compute some influence scores, and save the results in the `Influence` folder.
