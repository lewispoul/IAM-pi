from xtb_wrapper import run_xtb

xyz_path = "IAM_Results/methane_test.xyz"
output_file = "xtb_output/methane_xtb.log"

results = run_xtb(xyz_path, output_file)

for key, value in results.items():
    print(f"{key}: {value}")
