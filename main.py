import tkinter as tk
import tkinter.messagebox
from typing import Union, Callable
import customtkinter as ctk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import customtkinter
import math
import pandas as pd

customtkinter.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"

""" DataFrame Rectangular Plate for beta and alpha"""

columns_RP = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, np.inf]
data_RP = [[0.2874, 0.3762, 0.4530, 0.5172, 0.5688, 0.6102, 0.7134, 0.7410, 0.7476, 0.7500],
           [0.0444, 0.0616, 0.0770, 0.0906, 0.1017, 0.1110, 0.1335, 0.1400, 0.1417, 0.1421]]
index_RP = ['beta', 'alpha']
df_RP = pd.DataFrame(data_RP, columns=columns_RP, index=index_RP)

""" DataFrame Peak Velocity Pressure"""

columns_PVP = pd.MultiIndex.from_tuples(
    [('Area I', 'Coastal'), ('Area I', 'Rural'), ('Area I', 'Urban'),
     ('Area II', 'Coastal'), ('Area II', 'Rural'), ('Area II', 'Urban'),
     ('Area III', 'Rural'), ('Area III', 'Urban')],
    names=['Area', 'Type'])
index_PVP = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40,
             45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150]

data_PVP = np.array([
    [0.93, 0.71, 0.69, 0.78, 0.60, 0.58, 0.49, 0.48],
    [1.11, 0.71, 0.69, 0.93, 0.60, 0.58, 0.49, 0.48],
    [1.22, 0.71, 0.69, 1.02, 0.60, 0.58, 0.49, 0.48],
    [1.30, 0.71, 0.69, 1.09, 0.60, 0.58, 0.49, 0.48],
    [1.37, 0.78, 0.69, 1.14, 0.66, 0.58, 0.54, 0.48],
    [1.42, 0.84, 0.69, 1.19, 0.71, 0.58, 0.58, 0.48],
    [1.47, 0.89, 0.69, 1.23, 0.75, 0.58, 0.62, 0.48],
    [1.51, 0.94, 0.73, 1.26, 0.79, 0.62, 0.65, 0.51],
    [1.55, 0.98, 0.77, 1.29, 0.82, 0.65, 0.68, 0.53],
    [1.58, 1.02, 0.81, 1.32, 0.85, 0.68, 0.70, 0.56],
    [1.71, 1.16, 0.96, 1.43, 0.98, 0.80, 0.80, 0.66],
    [1.80, 1.27, 1.07, 1.51, 1.07, 0.90, 0.88, 0.74],
    [1.88, 1.36, 1.16, 1.57, 1.14, 0.97, 0.94, 0.80],
    [1.94, 1.43, 1.23, 1.63, 1.20, 1.03, 0.99, 0.85],
    [2.00, 1.50, 1.30, 1.67, 1.25, 1.09, 1.03, 0.89],
    [2.04, 1.55, 1.35, 1.71, 1.30, 1.13, 1.07, 0.93],
    [2.09, 1.60, 1.40, 1.75, 1.34, 1.17, 1.11, 0.97],
    [2.12, 1.65, 1.45, 1.78, 1.38, 1.21, 1.14, 1.00],
    [2.16, 1.69, 1.49, 1.81, 1.42, 1.25, 1.17, 1.03],
    [2.19, 1.73, 1.53, 1.83, 1.45, 1.28, 1.19, 1.05],
    [2.22, 1.76, 1.57, 1.86, 1.48, 1.31, 1.22, 1.08],
    [2.25, 1.80, 1.60, 1.88, 1.50, 1.34, 1.24, 1.10],
    [2.27, 1.83, 1.63, 1.90, 1.53, 1.37, 1.26, 1.13],
    [2.30, 1.86, 1.66, 1.92, 1.55, 1.39, 1.28, 1.15],
    [2.32, 1.88, 1.69, 1.94, 1.58, 1.42, 1.30, 1.17],
    [2.34, 1.91, 1.72, 1.96, 1.60, 1.44, 1.32, 1.18],
    [2.36, 1.93, 1.74, 1.98, 1.62, 1.46, 1.33, 1.20],
    [2.38, 1.96, 1.77, 1.99, 1.64, 1.48, 1.35, 1.22],
    [2.42, 2.00, 1.81, 2.03, 1.68, 1.52, 1.38, 1.25],
    [2.45, 2.04, 1.85, 2.05, 1.71, 1.55, 1.41, 1.28],
    [2.48, 2.08, 1.89, 2.08, 1.74, 1.59, 1.44, 1.31],
    [2.51, 2.12, 1.93, 2.10, 1.77, 1.62, 1.46, 1.33],
    [2.54, 2.15, 1.96, 2.13, 1.80, 1.65, 1.48, 1.35]])
df_PVP = pd.DataFrame(data=data_PVP, index=index_PVP, columns=columns_PVP)


class WidgetName(customtkinter.CTkFrame):
    def __init__(self, *args, width: int = 100, height: int = 32, **kwargs):
        super().__init__(*args, width=width, height=height, **kwargs)


class FloatSpinbox(customtkinter.CTkFrame):
    def __init__(self, *args, width: int = 100, height: int = 32,
                 step_size: Union[int, float] = 1, command: Callable = None, **kwargs):
        super().__init__(*args, width=width, height=height, **kwargs)

        self.step_size = step_size
        self.command = command

        self.configure(fg_color=("gray78", "gray28"))
        self.grid_columnconfigure((0, 2), weight=0)
        self.grid_columnconfigure(1, weight=1)

        self.subtract_button = customtkinter.CTkButton(
            self, text="-", width=height - 6, height=height - 6, command=self.subtract_button_callback)
        self.subtract_button.grid(row=0, column=0, padx=(3, 0), pady=3)

        self.entry = customtkinter.CTkEntry(self, width=width - (2 * height), height=height - 6, border_width=0)
        self.entry.grid(row=0, column=1, columnspan=1, padx=3, pady=3, sticky="ew")

        self.add_button = customtkinter.CTkButton(
            self, text="+", width=height - 6, height=height - 6, command=self.add_button_callback)
        self.add_button.grid(row=0, column=2, padx=(0, 3), pady=3)

        self.entry.insert(0, "0.0")

    def add_button_callback(self):
        if self.command is not None:
            self.command()
        try:
            value = float(self.entry.get()) + self.step_size
            self.entry.delete(0, "end")
            self.entry.insert(0, "{:.2f}".format(value))
        except ValueError:
            return

    def subtract_button_callback(self):
        if self.command is not None:
            self.command()
        try:
            value = float(self.entry.get()) - self.step_size
            self.entry.delete(0, "end")
            self.entry.insert(0, "{:.2f}".format(value))
        except ValueError:
            return

    def get(self) -> Union[float, None]:
        try:
            return float(self.entry.get())
        except ValueError:
            return None

    def set(self, value: float):
        self.entry.delete(0, "end")
        self.entry.insert(0, "{:.2f}".format(value))


class StructureHeightSpinBox(customtkinter.CTkFrame):
    def __init__(self, *args, width: int = 100, height: int = 32,
                 values: list = None, command: Callable = None, **kwargs):
        super().__init__(*args, width=width, height=height, **kwargs)

        self.values = values or []
        self.command = command

        self.configure(fg_color=("gray78", "gray28"))
        self.grid_columnconfigure((0, 1, 3, 4), weight=0)
        self.grid_columnconfigure(2, weight=1)

        self.jump_subtract_button = customtkinter.CTkButton(
            self, text="- -", width=height - 6, height=height - 6, command=self.jump_subtract_button_callback)
        self.jump_subtract_button.grid(row=0, column=0, padx=(3, 0), pady=3)

        self.subtract_button = customtkinter.CTkButton(
            self, text="-", width=height - 6, height=height - 6, command=self.subtract_button_callback)
        self.subtract_button.grid(row=0, column=1, padx=(3, 0), pady=3)

        self.entry_var = tk.StringVar(self)
        self.entry = customtkinter.CTkEntry(self, width=width - (2 * height), height=height - 6, border_width=0,
                                            state='disabled', textvariable=self.entry_var)
        self.entry.grid(row=0, column=2, columnspan=1, padx=3, pady=3, sticky="ew")

        self.add_button = customtkinter.CTkButton(
            self, text="+", width=height - 6, height=height - 6, command=self.add_button_callback)
        self.add_button.grid(row=0, column=3, padx=(0, 3), pady=3)

        self.jump_add_button = customtkinter.CTkButton(
            self, text="+ +", width=height - 6, height=height - 6, command=self.jump_add_button_callback)
        self.jump_add_button.grid(row=0, column=4, padx=(0, 3), pady=3)

        self.entry_var.set("50.0")

    def jump_subtract_button_callback(self):
        current_value = float(self.entry_var.get())
        if current_value in self.values:
            idx = self.values.index(current_value)
            if idx >= 5:
                self.entry_var.set("{:.1f}".format(self.values[idx - 5]))
                if self.command:
                    self.command()
        # self.update_button_states()

    def subtract_button_callback(self):
        current_value = float(self.entry_var.get())
        if current_value in self.values:
            idx = self.values.index(current_value)
            if idx > 0:
                self.entry_var.set("{:.1f}".format(self.values[idx - 1]))
                if self.command:
                    self.command()
        # self.update_button_states()

    def add_button_callback(self):
        current_value = float(self.entry_var.get())
        if current_value in self.values:
            idx = self.values.index(current_value)
            if idx < len(self.values) - 1:
                self.entry_var.set("{:.1f}".format(self.values[idx + 1]))
                if self.command:
                    self.command()
        # self.update_button_states()

    def jump_add_button_callback(self):
        current_value = float(self.entry_var.get())
        if current_value in self.values:
            idx = self.values.index(current_value)
            if idx < len(self.values) - 5:
                self.entry_var.set("{:.1f}".format(self.values[idx + 5]))
                if self.command:
                    self.command()
        # self.update_button_states()

    def get(self) -> Union[float, None]:
        try:
            return float(self.entry.get())
        except ValueError:
            return None

    


class DimensionWindow(ctk.CTkFrame):
    def __init__(self, master, df):
        super().__init__(master)

        self.df = df

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self.max_dimension_x = 1200
        self.max_dimension_y = 1500

        self.fig, self.ax = plt.subplots(figsize=(5, 5))
        self.ax.set_xlim(0, self.max_dimension_x)
        self.ax.set_ylim(0, self.max_dimension_y)
        self.ax.set_facecolor('lightslategray')
        self.ax.set_aspect('equal')

        self.rect = plt.Rectangle((0, 0), self.max_dimension_x / 2, self.max_dimension_y / 2, fc='lightsteelblue')
        self.ax.add_patch(self.rect)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # Sliders
        self.width_slider = tk.Scale(self, from_=0, to=self.max_dimension_x, length=400,
                                     command=self.update_rectangle, orient='horizontal')
        self.width_slider.set(self.max_dimension_x / 2)  # Set initial value to half of the maximum
        self.width_slider.grid(row=1, column=0, sticky="ew")

        self.height_slider = tk.Scale(self, from_=self.max_dimension_y, to=0, length=400,
                                      command=self.update_rectangle, orient='vertical')
        self.height_slider.set(self.max_dimension_y / 2)  # Set initial value to half of the maximum
        self.height_slider.grid(row=0, column=1, sticky="ns")

    def update_axes_limits(self):
        self.ax.set_xlim(0, self.max_dimension_x)
        self.ax.set_ylim(0, self.max_dimension_y)

    def update_rectangle(self, _):
        width = self.width_slider.get()
        height = self.height_slider.get()

        if width > height:
            width = height
            self.width_slider.set(width)

        self.rect.set_width(width)
        self.rect.set_height(height)
        self.fig.canvas.draw_idle()

        alpha, beta = self.calculate_ratio_and_get_alpha_beta(height, width)

    def calculate_ratio_and_get_alpha_beta(self, length, width):
        ratio = length / width if width != 0 else float('inf')

        # Find the closest column
        closest_column = min(self.df.columns,
                             key=lambda x: abs(x - ratio) if x != 'infinity' else abs(float('inf') - ratio))

        # Get alpha and beta
        alpha = self.df.loc['alpha', closest_column]
        beta = self.df.loc['beta', closest_column]

        return alpha, beta


class WindowParameters(customtkinter.CTkFrame):
    def __init__(self, master, dimension_window):
        super().__init__(master)

        self.master = master
        self.dimension_window = dimension_window

        self.identical_panes = tkinter.BooleanVar()

        # Variables for outer pane
        self.outer_strength_var = tkinter.StringVar()
        self.outer_thickness_var = tkinter.StringVar()
        self.outer_height_pane_var = tkinter.StringVar()
        self.outer_width_pane_var = tkinter.StringVar()

        # Outer pane
        self.frame_outer_pane = customtkinter.CTkFrame(self)
        self.frame_outer_pane.grid(row=1, column=0, sticky="nsew")

        self.label_outer_pane = customtkinter.CTkLabel(self.frame_outer_pane,
                                                       text="Outer Pane Window",
                                                       fg_color="gray28",
                                                       corner_radius=5)

        self.label_outer_pane.grid(row=0, column=0, columnspan=3)

        self.outer_strength_label = customtkinter.CTkLabel(self.frame_outer_pane, text="Strength:")
        self.outer_strength_label.grid(row=1, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.outer_strength_entry = customtkinter.CTkEntry(self.frame_outer_pane,
                                                           textvariable=self.outer_strength_var)
        self.outer_strength_entry.grid(row=1, column=1, sticky="w", padx=(5, 5), pady=(5, 5))
        self.outer_strength_unit = customtkinter.CTkLabel(self.frame_outer_pane, text="[N/mm^2]")
        self.outer_strength_unit.grid(row=1, column=2, sticky="w", padx=(5, 5), pady=(5, 5))

        self.outer_thickness_label = customtkinter.CTkLabel(self.frame_outer_pane, text="Thickness:")
        self.outer_thickness_label.grid(row=2, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.outer_thickness_entry = customtkinter.CTkEntry(self.frame_outer_pane,
                                                            textvariable=self.outer_thickness_var)
        self.outer_thickness_entry.grid(row=2, column=1, sticky="w", padx=(5, 5), pady=(5, 5))
        self.outer_thickness_unit = customtkinter.CTkLabel(self.frame_outer_pane, text="[mm]")
        self.outer_thickness_unit.grid(row=2, column=2, sticky="w", padx=(5, 5), pady=(5, 5))

        self.outer_height_pane_label = customtkinter.CTkLabel(self.frame_outer_pane, text="Height Pane (mm):")
        self.outer_height_pane_label.grid(row=3, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.outer_height_pane_entry = customtkinter.CTkEntry(self.frame_outer_pane,
                                                              textvariable=self.outer_height_pane_var)
        self.outer_height_pane_entry.grid(row=3, column=1, sticky="w", padx=(5, 5), pady=(5, 5))
        self.outer_height_pane_unit = customtkinter.CTkLabel(self.frame_outer_pane, text="[mm]")
        self.outer_height_pane_unit.grid(row=3, column=2, sticky="w", padx=(5, 5), pady=(5, 5))

        self.outer_width_pane_label = customtkinter.CTkLabel(self.frame_outer_pane, text="Width Pane (mm):")
        self.outer_width_pane_label.grid(row=4, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.outer_width_pane_entry = customtkinter.CTkEntry(self.frame_outer_pane,
                                                             textvariable=self.outer_width_pane_var)
        self.outer_width_pane_entry.grid(row=4, column=1, sticky="w", padx=(5, 5), pady=(5, 5))
        self.outer_width_pane_unit = customtkinter.CTkLabel(self.frame_outer_pane, text="[mm]")
        self.outer_width_pane_unit.grid(row=4, column=2, sticky="w", padx=(5, 5), pady=(5, 5))

        # Variables for inner pane
        self.inner_strength_var = tkinter.StringVar()
        self.inner_thickness_var = tkinter.StringVar()
        self.inner_height_pane_var = tkinter.StringVar()
        self.inner_width_pane_var = tkinter.StringVar()

        # Inner pane
        self.frame_inner_pane = customtkinter.CTkFrame(self)
        self.frame_inner_pane.grid(row=2, column=0, sticky="nsew")

        self.label_inner_pane = customtkinter.CTkLabel(self.frame_inner_pane,
                                                       text="Inner Pane Window",
                                                       fg_color="gray28",
                                                       corner_radius=5)
        self.label_inner_pane.grid(row=0, column=0, columnspan=3)

        self.inner_strength_label = customtkinter.CTkLabel(self.frame_inner_pane, text="Strength:")
        self.inner_strength_label.grid(row=1, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.inner_strength_entry = customtkinter.CTkEntry(self.frame_inner_pane,
                                                           textvariable=self.inner_strength_var)
        self.inner_strength_entry.grid(row=1, column=1, sticky="w", padx=(5, 5), pady=(5, 5))
        self.inner_strength_unit = customtkinter.CTkLabel(self.frame_inner_pane, text="[N/mm^2]")
        self.inner_strength_unit.grid(row=1, column=2, sticky="w", padx=(5, 5), pady=(5, 5))

        self.inner_thickness_label = customtkinter.CTkLabel(self.frame_inner_pane, text="Thickness:")
        self.inner_thickness_label.grid(row=2, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.inner_thickness_entry = customtkinter.CTkEntry(self.frame_inner_pane,
                                                            textvariable=self.inner_thickness_var)
        self.inner_thickness_entry.grid(row=2, column=1, sticky="w", padx=(5, 5), pady=(5, 5))
        self.inner_thickness_unit = customtkinter.CTkLabel(self.frame_inner_pane, text="[mm]")
        self.inner_thickness_unit.grid(row=2, column=2, sticky="w", padx=(5, 5), pady=(5, 5))

        self.inner_height_pane_label = customtkinter.CTkLabel(self.frame_inner_pane, text="Height Pane (mm):")
        self.inner_height_pane_label.grid(row=3, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.inner_height_pane_entry = customtkinter.CTkEntry(self.frame_inner_pane,
                                                              textvariable=self.inner_height_pane_var)
        self.inner_height_pane_entry.grid(row=3, column=1, sticky="w", padx=(5, 5), pady=(5, 5))
        self.inner_height_pane_unit = customtkinter.CTkLabel(self.frame_inner_pane, text="[mm]")
        self.inner_height_pane_unit.grid(row=3, column=2, sticky="w", padx=(5, 5), pady=(5, 5))

        self.inner_width_pane_label = customtkinter.CTkLabel(self.frame_inner_pane, text="Width Pane:")
        self.inner_width_pane_label.grid(row=4, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.inner_width_pane_entry = customtkinter.CTkEntry(self.frame_inner_pane,
                                                             textvariable=self.inner_width_pane_var)
        self.inner_width_pane_entry.grid(row=4, column=1, sticky="w", padx=(5, 5), pady=(5, 5))
        self.inner_width_pane_unit = customtkinter.CTkLabel(self.frame_inner_pane, text="[mm]")
        self.inner_width_pane_unit.grid(row=4, column=2, sticky="w", padx=(5, 5), pady=(5, 5))

        # Extra Parameters
        self.width_gap_var = tkinter.StringVar()
        self.k_a_var = tkinter.StringVar()
        self.k_a_var.set("1")  # Change this to the k_a function relying on width and height window
        self.k_e_var = tkinter.StringVar()
        self.k_e_var.set("1")
        self.k_mod_var = tkinter.StringVar()
        self.k_mod_var.set("1")
        self.k_sp_var = tkinter.StringVar()
        self.k_sp_var.set("1")
        self.gamma_var = tkinter.StringVar()
        self.gamma_var.set("1.6")

        self.frame_extra_params = customtkinter.CTkFrame(self)
        self.frame_extra_params.grid(row=3, column=0, sticky="nsew")

        self.label_extra_params = customtkinter.CTkLabel(self.frame_extra_params,
                                                         text="Extra Parameters",
                                                         fg_color="gray28",
                                                         corner_radius=5)
        self.label_extra_params.grid(row=0, column=0, columnspan=3)

        self.width_gap_label = customtkinter.CTkLabel(self.frame_extra_params, text="Width Gap:")
        self.width_gap_label.grid(row=1, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.width_gap_entry = customtkinter.CTkEntry(self.frame_extra_params,
                                                      textvariable=self.width_gap_var)
        self.width_gap_entry.grid(row=1, column=1, sticky="w", padx=(5, 5), pady=(5, 5))
        self.width_gap_unit = customtkinter.CTkLabel(self.frame_extra_params, text="[mm]")
        self.width_gap_unit.grid(row=1, column=2, sticky="w", padx=(5, 5), pady=(5, 5))

        self.k_a_label = customtkinter.CTkLabel(self.frame_extra_params, text="k_a:")
        self.k_a_label.grid(row=2, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.k_a_entry = customtkinter.CTkEntry(self.frame_extra_params,
                                                      textvariable=self.k_a_var)
        self.k_a_entry.grid(row=2, column=1, sticky="w", padx=(5, 5), pady=(5, 5))

        self.k_e_label = customtkinter.CTkLabel(self.frame_extra_params, text="k_e:")
        self.k_e_label.grid(row=3, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.k_e_entry = customtkinter.CTkEntry(self.frame_extra_params, textvariable=self.k_e_var)
        self.k_e_entry.grid(row=3, column=1, sticky="w", padx=(5, 5), pady=(5, 5))

        self.k_mod_label = customtkinter.CTkLabel(self.frame_extra_params, text="k_mod:")
        self.k_mod_label.grid(row=4, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.k_mod_entry = customtkinter.CTkEntry(self.frame_extra_params, textvariable=self.k_mod_var)
        self.k_mod_entry.grid(row=4, column=1, sticky="w", padx=(5, 5), pady=(5, 5))

        self.k_sp_label = customtkinter.CTkLabel(self.frame_extra_params, text="k_sp:")
        self.k_sp_label.grid(row=5, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.k_sp_entry = customtkinter.CTkEntry(self.frame_extra_params, textvariable=self.k_sp_var)
        self.k_sp_entry.grid(row=5, column=1, sticky="w", padx=(5, 5), pady=(5, 5))

        self.gamma_label = customtkinter.CTkLabel(self.frame_extra_params, text="gamma:")
        self.gamma_label.grid(row=6, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.gamma_entry = customtkinter.CTkEntry(self.frame_extra_params, textvariable=self.gamma_var)
        self.gamma_entry.grid(row=6, column=1, sticky="w", padx=(5, 5), pady=(5, 5))


        self.check_box = customtkinter.CTkCheckBox(self,
                                                   text="Inner and Outer Pane Identical",
                                                   variable=self.identical_panes,
                                                   command=self.update_pane_values)
        self.check_box.grid(row=0, column=0, sticky="nsew")

        # Callbacks
        self.outer_strength_var.trace('w', self.update_pane_values)
        self.outer_thickness_var.trace('w', self.update_pane_values)
        self.outer_height_pane_var.trace('w', self.update_pane_values)
        self.outer_width_pane_var.trace('w', self.update_pane_values)

        self.update_values_button = ctk.CTkButton(self, text="Update Values", command=self.update_values)
        self.update_values_button.grid(row=4, column=0)

    def update_pane_values(self, *args):
        if self.identical_panes.get():
            self.inner_strength_var.set(self.outer_strength_var.get())
            self.inner_thickness_var.set(self.outer_thickness_var.get())
            self.inner_height_pane_var.set(self.outer_height_pane_var.get())
            self.inner_width_pane_var.set(self.outer_width_pane_var.get())

            self.inner_strength_entry.configure(state='disabled')
            self.inner_thickness_entry.configure(state='disabled')
            self.inner_height_pane_entry.configure(state='disabled')
            self.inner_width_pane_entry.configure(state='disabled')
        else:
            self.inner_strength_entry.configure(state='normal')
            self.inner_thickness_entry.configure(state='normal')
            self.inner_height_pane_entry.configure(state='normal')
            self.inner_width_pane_entry.configure(state='normal')

    def update_values(self):
        self.master.values = {
            "outer_strength": self.outer_strength_var.get(),
            "outer_thickness": self.outer_thickness_var.get(),
            "outer_height_pane": self.outer_height_pane_var.get(),
            "outer_width_pane": self.outer_width_pane_var.get(),

            "inner_strength": self.inner_strength_var.get(),
            "inner_thickness": self.inner_thickness_var.get(),
            "inner_height_pane": self.inner_height_pane_var.get(),
            "inner_width_pane": self.inner_width_pane_var.get(),

            "width_gap": self.width_gap_var.get()
        }
        self.dimension_window.max_dimension_x = min(float(self.outer_width_pane_var.get()),
                                                    float(self.inner_width_pane_var.get()))
        self.dimension_window.max_dimension_y = min(float(self.outer_height_pane_var.get()),
                                                    float(self.inner_height_pane_var.get()))

        self.dimension_window.width_slider.configure(to=self.dimension_window.max_dimension_x)
        self.dimension_window.height_slider.configure(from_=self.dimension_window.max_dimension_y)

        self.dimension_window.update_axes_limits()
        print(self.master.values)

    def k_a(self):
        width = 1500 # change this to width of the rectangle on the left
        length = 2000 # change this to length of the rectangle on the left
        return 1.644 * (width * length) ** (-1.25)

class WindForce(customtkinter.CTkFrame):
    """" Class to set the parameters for the wind-pressure and calculate the wind-force"""

    def __init__(self, master):
        super().__init__(master)
        self.tab_data_dict = None
        self.wind_face_var = tk.StringVar()  # StringVar to keep track of wind-face selection

        self.tabview = customtkinter.CTkTabview(self, width=250)
        self.tabview.pack(fill='both', expand=True)
        self.tabview.add("Location")
        self.tabview.add("Wind-pressure")

        self.tabview.tab("Location").grid_columnconfigure(0, weight=1)
        self.tabview.tab("Wind-pressure").grid_columnconfigure(0, weight=1)

        self.optionmenu_1a = customtkinter.CTkOptionMenu(self.tabview.tab("Location"),
                                                         dynamic_resizing=False,
                                                         values=["AREA I", "AREA II", "AREA III"])
        self.optionmenu_1a.grid(row=0, column=0, padx=20, pady=(10, 10))

        self.optionmenu_1b = customtkinter.CTkOptionMenu(self.tabview.tab("Location"),
                                                         values=["Coastal", "Rural", "Urban"])
        self.optionmenu_1b.grid(row=1, column=0, padx=20, pady=(10, 10))

        self.label_height_1 = customtkinter.CTkLabel(self.tabview.tab("Location"),
                                                     text=f'Structure height [m]',
                                                     fg_color="gray28",
                                                     corner_radius=5)
        self.label_height_1.grid(row=2, column=0, padx=20, pady=(10, 0))
        list_height_structures = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,
                                  80,
                                  85, 90, 95, 100, 110, 120, 130, 140, 150]
        self.spinbox_1 = StructureHeightSpinBox(self.tabview.tab("Location"),
                                                width=150,
                                                values=list_height_structures)
        self.spinbox_1.grid(row=3, column=0, padx=20, pady=(0, 10))

        self.label_windface = customtkinter.CTkLabel(self.tabview.tab("Location"),
                                                     text=f'Wind-face',
                                                     fg_color="gray28",
                                                     corner_radius=5)
        self.label_windface.grid(row=4, column=0, padx=20, pady=(10, 0))

        self.optionmenu_1c = customtkinter.CTkOptionMenu(self.tabview.tab("Location"),
                                                         dynamic_resizing=False,
                                                         variable=self.wind_face_var,  # Use wind_face_var for value
                                                         values=["A", "B", "C", "D", "E"])
        self.optionmenu_1c.grid(row=5, column=0, padx=20, pady=(0, 10))

        self.button_1 = customtkinter.CTkButton(self.tabview.tab("Location"),
                                                text="Update values",
                                                command=self.update_values)
        self.button_1.grid(row=6, column=0, padx=20, pady=5)

        self.label_tab_2 = customtkinter.CTkLabel(self.tabview.tab("Wind-pressure"),
                                                  text="Wind-pressure []",
                                                  fg_color="gray28",
                                                  corner_radius=5)

        self.label_tab_2.grid(row=0, column=0, padx=20, pady=(20, 0))

        self.spinbox_2a = FloatSpinbox(self.tabview.tab("Wind-pressure"),
                                       width=150,
                                       step_size=0.1)
        self.spinbox_2a.grid(row=1, column=0, padx=20, pady=(0, 20))
        self.spinbox_2a.set(2.0)

        self.label_windface = customtkinter.CTkLabel(self.tabview.tab("Wind-pressure"),
                                                     text=f'Wind-face',
                                                     fg_color="gray28",
                                                     corner_radius=5)
        self.label_windface.grid(row=2, column=0, padx=20, pady=(10, 0))

        self.optionmenu_1c = customtkinter.CTkOptionMenu(self.tabview.tab("Wind-pressure"),
                                                         dynamic_resizing=False,
                                                         variable=self.wind_face_var,  # Use wind_face_var for value
                                                         values=["A", "B", "C", "D", "E"])
        self.optionmenu_1c.grid(row=3, column=0, padx=20, pady=(0, 10))

        self.button_2 = customtkinter.CTkButton(self.tabview.tab("Wind-pressure"),
                                                text="Update values",
                                                command=self.update_values)
        self.button_2.grid(row=4, column=0, padx=20, pady=5)

    def update_values(self):
        """ Update the values to the master library and switch out the old values"""

        selected_tab = self.tabview.get()
        if selected_tab == "Location":
            self.tab_data_dict = {
                "tab": selected_tab,
                "area": self.optionmenu_1a.get(),
                "location": self.optionmenu_1b.get(),
                "structure_height": self.spinbox_1.get(),
                "wind_face": self.wind_face_var.get()
            }
        elif selected_tab == "Wind-pressure":
            self.tab_data_dict = {
                "tab": selected_tab,
                "wind_pressure": self.spinbox_2a.get(),
                "wind_face": self.wind_face_var.get()
            }
        else:
            raise ValueError(f"Unexpected tab: {selected_tab}")
        print(self.tab_data_dict)


    # def


class ExternalLoadLending:
    """Class to compute and return the external load lending."""

    def __init__(self, side_length_short, side_length_long, gap_width, outer_thickness, inner_thickness, external_load):
        """Initialize the class with side lengths, gap width, thicknesses, and external load."""

        self.side_length_short = side_length_short
        self.side_length_long = side_length_long
        self.gap_width = gap_width
        self.outer_thickness = outer_thickness
        self.inner_thickness = inner_thickness
        self.external_load = external_load

    def compute_z1(self):
        """Compute and return the value of z1 based on the side lengths."""

        _exp = (-1.123 * ((self.side_length_long / self.side_length_short) - 1) ** 1.097)
        return 181.8 * (self.side_length_short / self.side_length_long) ** 2 * (
                0.00406 + 0.00896 * (1 - math.exp(_exp)))

    def compute_chi(self):
        """Compute and return the value of chi based on the side lengths and z1."""

        z1 = self.compute_z1()
        _exp = (-6.8 * (self.side_length_long / self.side_length_short) ** 1.33)
        return (z1 / 16) * (0.4198 + 0.22 * math.exp(_exp)) * (self.side_length_long / self.side_length_short) ** 2

    def compute_a_(self):
        """Compute and return the value of a_ based on the thicknesses, gap width, and chi."""

        chi = self.compute_chi()
        return 29.9 * ((self.gap_width * (self.outer_thickness ** 3) * (self.inner_thickness ** 3)) / (
                ((self.outer_thickness ** 3) + (self.inner_thickness ** 3)) * chi)) ** 0.25

    def compute_isolation_factor(self):
        """Compute and return the isolation factor based on the side lengths and a_."""

        a_ = self.compute_a_()
        return 1 / (1 + (self.side_length_short / a_) ** 4)

    def compute_pe(self):
        """Compute and return the value of P_E based on the isolation factor, thicknesses, and external load."""

        isolation_factor = self.compute_isolation_factor()
        return (1 - isolation_factor) * ((self.inner_thickness ** 3) / (
                (self.outer_thickness ** 3) + (self.inner_thickness ** 3))) * self.external_load


class UnityCheck(customtkinter.CTkFrame):
    def __init__(self, master):
        super().__init__(master)

        # Variables for outer pane

        self.UC_outer = customtkinter.CTkFrame(self)
        self.UC_outer.grid(row=1, column=0, sticky="nsew")

        self.label_outer_pane = customtkinter.CTkLabel(self.UC_outer,
                                                       text="U.C. outer pane",
                                                       fg_color="gray28",
                                                       corner_radius=5)

        self.label_outer_pane.grid(row=0, column=0, columnspan=2)

        self.ULS_outer_label = customtkinter.CTkLabel(self.UC_outer, text="ULS:")
        self.ULS_outer_label.grid(row=1, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.ULS_outer_entry = customtkinter.CTkEntry(self.UC_outer, width=80)
        self.ULS_outer_entry.grid(row=1, column=1, sticky="w", padx=(5, 5), pady=(5, 5))

        self.SLS_outer_label = customtkinter.CTkLabel(self.UC_outer, text="SLS:")
        self.SLS_outer_label.grid(row=2, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.SLS_outer_entry = customtkinter.CTkEntry(self.UC_outer, width=80)
        self.SLS_outer_entry.grid(row=2, column=1, sticky="w", padx=(5, 5), pady=(5, 5))

        self.UC_inner = customtkinter.CTkFrame(self)
        self.UC_inner.grid(row=2, column=0, sticky="nsew")

        self.label_inner_pane = customtkinter.CTkLabel(self.UC_inner,
                                                       text="U.C. inner pane",
                                                       fg_color="gray28",
                                                       corner_radius=5)

        self.label_inner_pane.grid(row=0, column=0, columnspan=2)

        self.ULS_inner_label = customtkinter.CTkLabel(self.UC_inner, text="ULS:")
        self.ULS_inner_label.grid(row=1, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.ULS_inner_entry = customtkinter.CTkEntry(self.UC_inner, width=80)
        self.ULS_inner_entry.grid(row=1, column=1, sticky="w", padx=(5, 5), pady=(5, 5))

        self.SLS_inner_label = customtkinter.CTkLabel(self.UC_inner, text="SLS:")
        self.SLS_inner_label.grid(row=2, column=0, sticky="w", padx=(5, 5), pady=(5, 5))
        self.SLS_inner_entry = customtkinter.CTkEntry(self.UC_inner, width=80)
        self.SLS_inner_entry.grid(row=2, column=1, sticky="w", padx=(5, 5), pady=(5, 5))

   
        

class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.title("Parametric tool")
        self.geometry(f"{1200}x{680}")

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure((1, 2, 3), weight=0)
        self.grid_rowconfigure((0, 1, 2), weight=1)

        self.dimensions_window_frame = DimensionWindow(self, df_RP)
        self.dimensions_window_frame.grid(row=0, rowspan=2, column=0, columnspan=2, padx=(10, 0), pady=(20, 0),
                                          sticky="n")

        self.windforce_frame = WindForce(self)
        self.windforce_frame.grid(row=0, column=4, padx=(10, 0), pady=(0, 0), sticky="n")

        self.window_parameters_frame = WindowParameters(self, self.dimensions_window_frame)
        self.window_parameters_frame.grid(row=0, rowspan=2, column=3, padx=(10, 0), pady=0, sticky="nsew")

        self.unity_check_frame = UnityCheck(self)
        self.unity_check_frame.grid(row=1, column=4, padx=(10, 0), pady=(10, 0), ipadx=90, sticky='s')


if __name__ == "__main__":
    app = App()
    app.mainloop()