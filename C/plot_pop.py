import pandas as pd
import matplotlib.pyplot as plt
import os
#print(os.getcwd())
# Read the CSV file
#
data = "c1a"
data_folder = os.path.join(os.path.dirname(__file__), "Solutions", data)
df = pd.read_csv(os.path.join(data_folder, f"{data}_population.csv"))

# Plot populations vs time
plt.figure()
plt.plot(df["time"], df["population_u"], label="population_u")
plt.plot(df["time"], df["population_v"], label="population_v")
plt.plot(df["time"], df["population_w"], label="population_w")

# Labels and legend
plt.xlabel("Time")
plt.ylabel("Population")
plt.title("Populations vs Time, u0 = 0")
plt.legend()
plt.grid(True)

#plt.show()
plt.savefig(os.path.join(os.path.dirname(__file__), "Figures", f"{data}_populations.png"), dpi=150)