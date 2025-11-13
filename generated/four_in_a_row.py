import sys

ROWS = 6
COLS = 7

def create_board():
    """
    Creates and returns an empty game board as a list of lists.
    0 represents an empty slot, 1 for Player 1's piece, 2 for Player 2's piece.
    The board is structured such that board[0][0] is the bottom-left slot.
    """
    return [[0 for _ in range(COLS)] for _ in range(ROWS)]

def print_board(board):
    """
    Prints the current state of the game board to the console.
    'X' represents Player 1, 'O' represents Player 2.
    """
    print("\n")
    # Print from top row down to simulate gravity
    for r in range(ROWS - 1, -1, -1):
        row_str = "|"
        for c in range(COLS):
            if board[r][c] == 0:
                row_str += "   |"
            elif board[r][c] == 1:
                row_str += " X |"
            else: # board[r][c] == 2
                row_str += " O |"
        print(row_str)
    print("----------------------------")
    # Print column numbers for player input
    col_numbers = "  "
    for i in range(1, COLS + 1):
        col_numbers += f"{i}   "
    print(col_numbers.strip())
    sys.stdout.flush() # Ensure output is immediately visible

def is_valid_location(board, col):
    """
    Checks if the given column is a valid location to drop a piece.
    A location is valid if the top-most row of that column is not yet filled.
    """
    return 0 <= col < COLS and board[ROWS - 1][col] == 0

def get_next_open_row(board, col):
    """
    Finds the lowest empty row in the given column where a piece can be dropped.
    """
    for r in range(ROWS):
        if board[r][col] == 0:
            return r
    return -1 # Should not happen if is_valid_location is checked first

def drop_piece(board, row, col, piece):
    """
    Places the given player's piece (1 or 2) into the specified row and column on the board.
    """
    board[row][col] = piece

def winning_move(board, piece):
    """
    Checks if the given player (piece) has made a winning move.
    Checks for 4 consecutive pieces horizontally, vertically, and diagonally.
    """
    # Check horizontal locations for win
    for c in range(COLS - 3):
        for r in range(ROWS):
            if all(board[r][c+i] == piece for i in range(4)):
                return True

    # Check vertical locations for win
    for c in range(COLS):
        for r in range(ROWS - 3):
            if all(board[r+i][c] == piece for i in range(4)):
                return True

    # Check positively sloped diagonals (bottom-left to top-right)
    for c in range(COLS - 3):
        for r in range(ROWS - 3):
            if all(board[r+i][c+i] == piece for i in range(4)):
                return True

    # Check negatively sloped diagonals (top-left to bottom-right)
    for c in range(COLS - 3):
        for r in range(3, ROWS): # Start from row 3 (0-indexed) to go up
            if all(board[r-i][c+i] == piece for i in range(4)):
                return True
    return False

def is_board_full(board):
    """
    Checks if the game board is completely filled (resulting in a tie if no winner).
    """
    return all(board[ROWS - 1][c] != 0 for c in range(COLS))

def four_in_a_row_game():
    """
    Plays a command-line based Four-in-a-Row game between two players.
    """
    board = create_board()
    game_over = False
    turn = 0 # 0 for Player 1, 1 for Player 2

    print("Welcome to Four-in-a-Row!")
    print_board(board)

    while not game_over:
        if turn == 0:
            player_piece = 1
            player_name = "Player 1 (X)"
        else:
            player_piece = 2
            player_name = "Player 2 (O)"

        while True:
            try:
                col_input = input(f"{player_name}, choose a column (1-{COLS}): ")
                col = int(col_input) - 1 # Convert to 0-indexed column
                
                if not (0 <= col < COLS):
                    print(f"Column must be between 1 and {COLS}. Please try again.")
                elif not is_valid_location(board, col):
                    print(f"Column {col + 1} is full. Please choose another column.")
                else:
                    break
            except ValueError:
                print("Invalid input. Please enter a number.")
            except EOFError: # Handle Ctrl+D or unexpected end of input
                print("\nGame aborted by user.")
                return

        row = get_next_open_row(board, col)
        drop_piece(board, row, col, player_piece)

        print_board(board)

        if winning_move(board, player_piece):
            print(f"\n{player_name} wins!")
            game_over = True
        elif is_board_full(board):
            print("\nIt's a tie!")
            game_over = True

        turn = (turn + 1) % 2 # Switch turns

    print("Game Over!")