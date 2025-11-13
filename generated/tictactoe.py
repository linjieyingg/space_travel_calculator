import random

def play_tic_tac_toe():
    """
    Plays a game of Tic-Tac-Toe between two players (X and O).
    The game is played on a 3x3 board. Players take turns marking
    empty spaces. The first player to get three of their marks
    in a row (horizontal, vertical, or diagonal) wins.
    If all spaces are filled and no player has won, the game is a draw.
    """

    board = [' ' for _ in range(9)]
    player_symbols = ['X', 'O']
    current_player_idx = random.randint(0, 1) # Randomly choose starting player
    current_player_symbol = player_symbols[current_player_idx]

    def display_board(current_board):
        """Prints the Tic-Tac-Toe board to the console."""
        print(f"\n {current_board[0]} | {current_board[1]} | {current_board[2]} ")
        print("---+---+---")
        print(f" {current_board[3]} | {current_board[4]} | {current_board[5]} ")
        print("---+---+---")
        print(f" {current_board[6]} | {current_board[7]} | {current_board[8]} \n")

    def check_win(current_board, player):
        """Checks if the given player has won the game."""
        win_conditions = [
            # Horizontal
            [0, 1, 2], [3, 4, 5], [6, 7, 8],
            # Vertical
            [0, 3, 6], [1, 4, 7], [2, 5, 8],
            # Diagonal
            [0, 4, 8], [2, 4, 6]
        ]
        for condition in win_conditions:
            if all(current_board[i] == player for i in condition):
                return True
        return False

    def check_draw(current_board):
        """Checks if the board is full (a draw)."""
        return ' ' not in current_board

    print("Welcome to Tic-Tac-Toe!")
    print("Players will choose 'X' or 'O'.")
    print("Enter a number from 1-9 to place your mark on the board:")
    display_board(['1', '2', '3', '4', '5', '6', '7', '8', '9'])
    print(f"Player {current_player_symbol} goes first!\n")

    game_over = False
    while not game_over:
        display_board(board)
        print(f"It's Player {current_player_symbol}'s turn.")

        valid_move = False
        while not valid_move:
            try:
                move = int(input("Enter your move (1-9): ")) - 1
                if 0 <= move <= 8 and board[move] == ' ':
                    board[move] = current_player_symbol
                    valid_move = True
                else:
                    print("Invalid move. That spot is either taken or out of range (1-9). Please try again.")
            except ValueError:
                print("Invalid input. Please enter a number between 1 and 9.")

        if check_win(board, current_player_symbol):
            display_board(board)
            print(f"Congratulations! Player {current_player_symbol} wins!")
            game_over = True
        elif check_draw(board):
            display_board(board)
            print("It's a draw!")
            game_over = True
        else:
            # Switch player
            current_player_idx = (current_player_idx + 1) % 2
            current_player_symbol = player_symbols[current_player_idx]

    print("Game Over!")

if __name__ == '__main__':
    play_tic_tac_toe()