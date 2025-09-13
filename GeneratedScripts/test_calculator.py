import unittest
from unittest.mock import MagicMock, patch
import tkinter as tk
from calculator import Calculator

class TestCalculator(unittest.TestCase):
    def setUp(self):
        # Set up a Tkinter root window
        self.root = tk.Tk()
        self.calculator = Calculator(self.root)

    def tearDown(self):
        # Destroy the Tkinter root window after each test
        self.root.destroy()

    def test_initial_state(self):
        # Test initial state of the calculator
        self.assertEqual(self.calculator.expression, "")
        self.assertEqual(self.calculator.input_text.get(), "")

    def test_button_click_number(self):
        # Simulate clicking number buttons
        self.calculator.on_button_click('1')
        self.calculator.on_button_click('2')
        self.calculator.on_button_click('3')
        self.assertEqual(self.calculator.expression, "123")
        self.assertEqual(self.calculator.input_text.get(), "123")

    def test_button_click_operator(self):
        # Simulate clicking operator buttons
        self.calculator.on_button_click('1')
        self.calculator.on_button_click('+')
        self.calculator.on_button_click('2')
        self.assertEqual(self.calculator.expression, "1+2")
        self.assertEqual(self.calculator.input_text.get(), "1+2")

    def test_calculate_result(self):
        # Test calculating a valid expression
        self.calculator.expression = "2+3*4"
        self.calculator.calculate_result()
        self.assertEqual(self.calculator.expression, "14")
        self.assertEqual(self.calculator.input_text.get(), "14")

    def test_calculate_result_invalid(self):
        # Test calculating an invalid expression
        self.calculator.expression = "2+*3"
        with patch('tkinter.messagebox.showerror') as mock_showerror:
            self.calculator.calculate_result()
            mock_showerror.assert_called_once_with("Error", "Invalid Input")
        self.assertEqual(self.calculator.expression, "")
        self.assertEqual(self.calculator.input_text.get(), "")

    def test_clear_input(self):
        # Test clearing the input
        self.calculator.expression = "123"
        self.calculator.input_text.set("123")
        self.calculator.clear_input()
        self.assertEqual(self.calculator.expression, "")
        self.assertEqual(self.calculator.input_text.get(), "")

if __name__ == '__main__':
    unittest.main()