import numpy as np
import tkinter as tk
from tkinter import messagebox
from scipy.optimize import linprog

class TransportationOptimizer:
    def __init__(self, root):  # &lt;-- Corrected: double underscores
        self.root = root
        self.root.title("Transportation Problem Optimizer")
        self.root.geometry("500x400")
        
        self.create_widgets()

    def create_widgets(self):
        tk.Label(self.root, text="Number of Sources:").grid(row=0, column=0)
        self.sources_entry = tk.Entry(self.root)
        self.sources_entry.grid(row=0, column=1)
        
        tk.Label(self.root, text="Number of Destinations:").grid(row=1, column=0)
        self.destinations_entry = tk.Entry(self.root)
        self.destinations_entry.grid(row=1, column=1)
        
        tk.Button(self.root, text="Enter Data", command=self.get_inputs).grid(row=2, column=0, columnspan=2)

    def get_inputs(self):
        try:
            self.num_sources = int(self.sources_entry.get())
            self.num_destinations = int(self.destinations_entry.get())
            self.input_window()
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers")

    def input_window(self):
        self.input_win = tk.Toplevel(self.root)
        self.input_win.title("Enter Cost, Supply, and Demand")
        
        self.cost_entries = []
        self.supply_entries = []
        self.demand_entries = []
        
        tk.Label(self.input_win, text="Cost Matrix").grid(row=0, column=1, columnspan=self.num_destinations)
        for i in range(self.num_sources):
            row_entries = []
            for j in range(self.num_destinations):
                entry = tk.Entry(self.input_win, width=5)
                entry.grid(row=i+1, column=j+1)
                row_entries.append(entry)
            self.cost_entries.append(row_entries)
        
        tk.Label(self.input_win, text="Supply").grid(row=0, column=self.num_destinations+1)
        for i in range(self.num_sources):
            entry = tk.Entry(self.input_win, width=5)
            entry.grid(row=i+1, column=self.num_destinations+1)
            self.supply_entries.append(entry)
        
        tk.Label(self.input_win, text="Demand").grid(row=self.num_sources+1, column=1, columnspan=self.num_destinations)
        for j in range(self.num_destinations):
            entry = tk.Entry(self.input_win, width=5)
            entry.grid(row=self.num_sources+2, column=j+1)
            self.demand_entries.append(entry)
        
        tk.Button(self.input_win, text="Solve", command=self.solve_transportation).grid(row=self.num_sources+3, column=0, columnspan=self.num_destinations+2)

    def solve_transportation(self):
        try:
            cost_matrix = np.array([[float(entry.get()) for entry in row] for row in self.cost_entries])
            supply = np.array([float(entry.get()) for entry in self.supply_entries])
            demand = np.array([float(entry.get()) for entry in self.demand_entries])
            
            num_supply, num_demand = cost_matrix.shape
            
            c = cost_matrix.flatten()
            A_eq = []
            b_eq = []
            
            for i in range(num_supply):
                row = [0] * (num_supply * num_demand)
                for j in range(num_demand):
                    row[i * num_demand + j] = 1
                A_eq.append(row)
                b_eq.append(supply[i])
            
            for j in range(num_demand):
                row = [0] * (num_supply * num_demand)
                for i in range(num_supply):
                    row[i * num_demand + j] = 1
                A_eq.append(row)
                b_eq.append(demand[j])
            
            result = linprog(c, A_eq=A_eq, b_eq=b_eq, method='highs')
            if result.success:
                allocations = result.x.reshape((num_supply, num_demand))
                allocation_str = "\n".join(["\t".join([f"{val:.2f}" for val in row]) for row in allocations])
                messagebox.showinfo("Solution", f"Optimal Cost: {result.fun:.2f}\n\nAllocations:\n{allocation_str}")
            else:
                messagebox.showerror("Error", "No solution found")
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numerical values")

if __name__ == "__main__":  
    root = tk.Tk()
    app = TransportationOptimizer(root)
    root.mainloop()