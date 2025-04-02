import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file (replace with your filename)
df = pd.read_csv("../generations/scores.csv")

# Convert generation and score to numeric (if not already)
df['generation'] = pd.to_numeric(df['generation'], errors='coerce')
df['score'] = pd.to_numeric(df['score'], errors='coerce')

# Group by generation and calculate best score per generation
best_scores = df.groupby('generation')['score'].min()

# Plot
plt.figure(figsize=(10, 6))
plt.plot(best_scores.index, best_scores.values, marker='o')
plt.title('Best Score per Generation')
plt.xlabel('Generation')
plt.ylabel('Score')
plt.grid(True)
plt.tight_layout()
plt.show()
