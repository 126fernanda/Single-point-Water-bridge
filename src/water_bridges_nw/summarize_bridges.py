import json
import pandas as pd
import sys

def summarize_jsonl(input_file, output_file):
    data = []
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            record = json.loads(line)
            
            if record.get('type') == 'frame':
                frame_idx = record['frame_idx']
                for path in record.get('paths', []):
                    data.append({
                        'Frame': frame_idx,
                        'Length': path['length'],
                        'Probability': round(path['probability'], 4),
                        'Avg_OO_Dist (Å)': round(path['avg_oo_dist'], 3),
                        'Atom_IDs': " -> ".join(map(str, path['atom_ids']))
                    })
                    
    if not data:
        print("No pathways found in the file.")
        return

    df = pd.DataFrame(data)
    df = df.sort_values(by=['Frame', 'Probability'], ascending=[True, False])
    
    with open(output_file, 'w') as out:
        # 1. GLOBAL SUMMARY
        out.write("=================== GLOBAL SUMMARY ===================\n")
        summary = df['Length'].value_counts().sort_index()
        for length, count in summary.items():
            out.write(f"Total water-bridges of length {length}: {count}\n")
        out.write("======================================================\n\n")

        # 2. PER-FRAME SUMMARY (Matrix format)
        out.write("================== PER-FRAME SUMMARY =================\n")
        out.write("Reads as: Number of bridges of each length per frame\n\n")
        
        # Create a pivot table: Rows=Frame, Columns=Length, Values=Count
        frame_summary = df.groupby(['Frame', 'Length']).size().unstack(fill_value=0)
        frame_summary.columns = [f"Len_{col}" for col in frame_summary.columns]
        frame_summary['Total_Bridges'] = frame_summary.sum(axis=1)
        
        out.write(frame_summary.to_string() + "\n")
        out.write("======================================================\n\n")

        # 3. DETAILED PATHWAYS
        out.write("====================== DETAILED PATHWAYS ======================\n")
        out.write(df.to_string(index=False) + "\n")
        out.write("===============================================================\n")
        
    print(f"Summary successfully saved to {output_file}")

if __name__ == "__main__":
    input_f = sys.argv[1] if len(sys.argv) > 1 else "results.jsonl"
    output_f = sys.argv[2] if len(sys.argv) > 2 else "summarize-water-bridges.txt"
    summarize_jsonl(input_f, output_f)