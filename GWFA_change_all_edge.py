import os
from GWFA_truncate import truncate 

NUM_EDGES = 12
BASE_DIR = "."  

def is_target_dir(name: str) -> bool:
    return name.startswith("out_trim_err") 

def main():
    for root, dirs, files in os.walk(BASE_DIR):
        if is_target_dir(os.path.basename(root)):
            for fn in files:
                if fn.endswith(".gfa"):
                    gfa_path = os.path.join(root, fn)
                    try:
                        truncate(gfa_path, NUM_EDGES)
                    except Exception as e:
                        print(f"[FAIL] {gfa_path}: {e}")

if __name__ == "__main__":
    main()
