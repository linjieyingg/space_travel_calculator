import numpy as np

def four_in_a_row():
    """
    Implements a console-based 4-in-a-row game.
    Players take turns dropping their pieces into a column.
    The first player to get four of their pieces in a row (horizontally,
    vertically, or diagonally) wins.
    """

    ROW_COUNT = 6
    COLUMN_COUNT = 7
    PLAYER_PIECES = {1: 'X', 2: 'O'} # Player 1 uses 'X', Player 2 uses 'O'
    EMPTY_SLOT = ' '

    def create_board():
        """Creates an empty 4-in-a-row game board."""
        return np.full((ROW_COUNT, COLUMN_COUNT), EMPTY_SLOT)

    def print_board(board):
        """Prints the game board to the console."""
        print("\n" + "-" * (COLUMN_COUNT * 4 + 1))
        for r in range(ROW_COUNT - 1, -1, -1): # Print from top to bottom logically
            row_str = " | ".join(board[r, :])
            print(f"| {row_str} |")
            print("-" * (COLUMN_COUNT * 4 + 1))
        print("   " + "   ".join(map(str, range(COLUMN_COUNT)))) # Column numbers
        print("\n")

    def is_valid_location(board, col):
        """Checks if a column is a valid place to drop a piece."""
        return 0 <= col < COLUMN_COUNT and board[ROW_COUNT - 1, col] == EMPTY_SLOT

    def get_next_open_row(board, col):
        """Finds the lowest open row in a given column."""
        for r in range(ROW_COUNT):
            if board[r, col] == EMPTY_SLOT:
                return r
        return -1 # Should not happen if is_valid_location is checked first

    def drop_piece(board, row, col, piece):
        """Places a piece on the board."""
        board[row, col] = piece

    def check_win(board, piece):
        """Checks if the given player has won the game."""

        # Check horizontal locations for win
        for c in range(COLUMN_COUNT - 3):
            for r in range(ROW_COUNT):
                if all(board[r, c+i] == piece for i in range(4)):
                    return True

        # Check vertical locations for win
        for c in range(COLUMN_COUNT):
            for r in range(ROW_COUNT - 3):
                if all(board[r+i, c] == piece for i in range(4)):
                    return True

        # Check positively sloped diagonals
        for c in range(COLUMN_COUNT - 3):
            for r in range(ROW_COUNT - 3):
                if all(board[r+i, c+i] == piece for i in range(4)):
                    return True

        # Check negatively sloped diagonals
        for c in range(COLUMN_COUNT - 3):
            for r in range(3, ROW_COUNT):
                if all(board[r-i, c+i] == piece for i in range(4)):
                    return True
        return False

    def is_board_full(board):
        """Checks if the board is full (for a draw)."""
        return EMPTY_SLOT not in board[ROW_COUNT - 1, :]

    board = create_board()
    game_over = False
    turn = 0 # 0 for Player 1, 1 for Player 2

    print("Welcome to Four in a Row!")
    print_board(board)

    while not game_over:
        player = (turn % 2) + 1
        current_piece = PLAYER_PIECES[player]

        try:
            col_input = input(f"Player {player} ({current_piece}), choose a column (0-{COLUMN_COUNT-1}): ")
            col = int(col_input)
        except ValueError:
            print("Invalid input. Please enter a number.")
            continue

        if is_valid_location(board, col):
            row = get_next_open_row(board, col)
            drop_piece(board, row, col, current_piece)

            print_board(board)

            if check_win(board, current_piece):
                print(f"Player {player} ({current_piece}) wins! Congratulations!")
                game_over = True
            elif is_board_full(board):
                print("The board is full! It's a draw!")
                game_over = True
            else:
                turn += 1
        else:
            print(f"Column {col} is invalid or full. Please choose another column.")

if __name__ == '__main__':
    four_in_a_row()