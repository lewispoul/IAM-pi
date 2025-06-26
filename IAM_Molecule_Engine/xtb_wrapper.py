
import subprocess

def run_xtb(input_file, output_file="xtb_output.log"):
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
    except subprocess.CalledProcessError as e:
        print("❌ XTB failed:")
        print(e.stderr)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run XTB on an XYZ file.")
    parser.add_argument("input", help="Path to the input XYZ file")
    parser.add_argument("--output", default="xtb_output.log", help="Output log file")
    args = parser.parse_args()

    run_xtb(args.input, args.output)
