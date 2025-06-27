import subprocess
import os

def run_xtb(input_file, output_file="xtb_output/xtb_output.log"):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    try:
        result = subprocess.run(
            ["xtb", input_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        with open(output_file, "w") as f:
            f.write(result.stdout)
        print(f"✅ XTB calculation complete. Output saved to: {output_file}")
        return parse_xtb_output(output_file)
    except subprocess.CalledProcessError as e:
        print("❌ XTB execution failed:")
        print(e.stderr)
        return {}

def parse_xtb_output(output_file):
    results = {}
    with open(output_file, "r") as f:
        for line in f:
            if "TOTAL ENERGY" in line:
                try:
                    results["energy"] = float(line.split()[-2])
                except:
                    results["energy"] = "⚠️ parse error"
            elif "HOMO-LUMO GAP" in line:
                try:
                    results["homo_lumo_gap"] = float(line.split()[-2])
                except:
                    results["homo_lumo_gap"] = "⚠️ parse error"
            elif "molecular dipole:" in line:
                next_line = next(f)
                parts = next_line.strip().split()
                if len(parts) >= 4:
                    results["dipole_x"] = float(parts[0])
                    results["dipole_y"] = float(parts[1])
                    results["dipole_z"] = float(parts[2])
                    results["dipole_total"] = float(parts[3])
    return results

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run XTB on an XYZ file.")
    parser.add_argument("input", help="Path to the input XYZ file")
    parser.add_argument("--output", default="xtb_output/xtb_output.log", help="Output log file")
    args = parser.parse_args()

    results = run_xtb(args.input, args.output)
    print(results)
