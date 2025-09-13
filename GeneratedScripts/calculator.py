import tkinter as tk
from tkinter import messagebox

class Calculator:
    def __init__(self, root):
        """
        Initialize the calculator with a Tkinter root window.
        
        :param root: The root window for the calculator GUI.
        """
        self.root = root
        self.root.title("Calculator")
        self.root.geometry("400x600")
        self.root.resizable(0, 0)
        
        self.expression = ""
        self.input_text = tk.StringVar()
        
        self.create_widgets()

    def create_widgets(self):
        """
        Create the widgets for the calculator GUI.
        """
        input_frame = tk.Frame(self.root, width=400, height=50, bd=0, highlightbackground="black", highlightcolor="black", highlightthickness=1)
        input_frame.pack(side=tk.TOP)
        
        input_field = tk.Entry(input_frame, font=('arial', 18, 'bold'), textvariable=self.input_text, width=50, bg="#eee", bd=0, justify=tk.RIGHT)
        input_field.grid(row=0, column=0)
        input_field.pack(ipady=10)
        
        btns_frame = tk.Frame(self.root, width=400, height=450, bg="grey")
        btns_frame.pack()
        
        self.create_buttons(btns_frame)

    def create_buttons(self, frame):
        """
        Create the buttons for the calculator and attach them to the frame.
        
        :param frame: The frame where the buttons will be placed.
        """
        buttons = [
            ('7', 1, 0), ('8', 1, 1), ('9', 1, 2), ('/', 1, 3),
            ('4', 2, 0), ('5', 2, 1), ('6', 2, 2), ('*', 2, 3),
            ('1', 3, 0), ('2', 3, 1), ('3', 3, 2), ('-', 3, 3),
            ('C', 4, 0), ('0', 4, 1), ('=', 4, 2), ('+', 4, 3),
        ]
        
        for (text, row, col) in buttons:
            button = tk.Button(frame, text=text, fg="black", width=10, height=3, bd=0, bg="#fff", cursor="hand2",
                               command=lambda t=text: self.on_button_click(t))
            button.grid(row=row, column=col, padx=1, pady=1)

    def on_button_click(self, char):
        """
        Handle button click events.
        
        :param char: The character of the button that was clicked.
        """
        if char == "=":
            self.calculate_result()
        elif char == "C":
            self.clear_input()
        else:
            self.expression += str(char)
            self.input_text.set(self.expression)

    def calculate_result(self):
        """
        Calculate the result of the current expression.
        """
        try:
            result = str(eval(self.expression))
            self.input_text.set(result)
            self.expression = result
        except Exception as e:
            messagebox.showerror("Error", "Invalid Input")
            self.expression = ""
            self.input_text.set("")

    def clear_input(self):
        """
        Clear the current input.
        """
        self.expression = ""
        self.input_text.set("")

def main():
    """
    The main function to run the calculator application.
    """
    root = tk.Tk()
    Calculator(root)
    root.mainloop()

if __name__ == "__main__":
    main()