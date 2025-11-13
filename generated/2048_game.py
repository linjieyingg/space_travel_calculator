import random
import os

class Game2048:
    def __init__(self):
        self.grid = [[0] * 4 for _ in range(4)]
        self.score = 0
        self.add_random_tile()
        self.add_random_tile()
        self.won = False
        self.game_over = False

    def add_random_tile(self):
        empty_cells = [(r, c) for r in range(4) for c in range(4) if self.grid[r][c] == 0]
        if empty_cells:
            r, c = random.choice(empty_cells)
            self.grid[r][c] = 2 if random.random() < 0.9 else 4 # 90% chance of 2, 10% chance of 4
            return True
        return False

    def display_board(self):
        # Clear console (works on Windows and Unix-like systems)
        os.system('cls' if os.name == 'nt' else 'clear') 
        print(f"Score: {self.score}")
        print("-" * 25)
        for row in self.grid:
            print("|", end="")
            for cell in row:
                if cell == 0:
                    print(f"{' ':^5}|", end="")
                else:
                    print(f"{cell:^5}|", end="")
            print()
            print("-" * 25)

    def _get_row_view(self, r):
        return self.grid[r]

    def _set_row_view(self, r, new_row):
        self.grid[r] = new_row

    def _get_col_view(self, c):
        return [self.grid[r][c] for r in range(4)]

    def _set_col_view(self, c, new_col):
        for r in range(4):
            self.grid[r][c] = new_col[r]

    def _shift_and_merge(self, line):
        original_line = list(line) # Make a copy to compare later

        # 1. Shift non-zero numbers to the left
        new_line = [num for num in line if num != 0]
        missing_zeros = [0] * (4 - len(new_line))
        new_line.extend(missing_zeros)

        # 2. Merge adjacent identical numbers
        for i in range(3):
            if new_line[i] != 0 and new_line[i] == new_line[i+1]:
                new_line[i] *= 2
                self.score += new_line[i]
                if new_line[i] == 2048:
                    self.won = True
                new_line[i+1] = 0 # Set merged tile to 0

        # 3. Shift again after merging
        final_line = [num for num in new_line if num != 0]
        missing_zeros_final = [0] * (4 - len(final_line))
        final_line.extend(missing_zeros_final)

        return final_line, final_line != original_line # Return new line and if it changed

    def move(self, direction):
        moved_happened = False

        if direction == 'up':
            for c in range(4):
                col = self._get_col_view(c)
                new_col, changed = self._shift_and_merge(col)
                if changed:
                    self._set_col_view(c, new_col)
                    moved_happened = True
        elif direction == 'down':
            for c in range(4):
                col = self._get_col_view(c)[::-1] # Reverse for down
                new_col, changed = self._shift_and_merge(col)
                if changed:
                    self._set_col_view(c, new_col[::-1]) # Reverse back
                    moved_happened = True
        elif direction == 'left':
            for r in range(4):
                row = self._get_row_view(r)
                new_row, changed = self._shift_and_merge(row)
                if changed:
                    self._set_row_view(r, new_row)
                    moved_happened = True
        elif direction == 'right':
            for r in range(4):
                row = self._get_row_view(r)[::-1] # Reverse for right
                new_row, changed = self._shift_and_merge(row)
                if changed:
                    self._set_row_view(r, new_row[::-1]) # Reverse back
                    moved_happened = True
        else:
            # Invalid direction, no change in game state
            return False

        if moved_happened:
            self.add_random_tile()
            self.check_game_over()
        return moved_happened

    def check_game_over(self):
        # If there are empty cells, game is not over
        if any(0 in row for row in self.grid):
            self.game_over = False
            return False

        # Check for possible merges in rows
        for r in range(4):
            for c in range(3):
                if self.grid[r][c] == self.grid[r][c+1]:
                    self.game_over = False
                    return False

        # Check for possible merges in columns
        for c in range(4):
            for r in range(3):
                if self.grid[r][c] == self.grid[r+1][c]:
                    self.game_over = False
                    return False

        # No empty cells and no possible merges, game is over
        self.game_over = True
        return True

def play_2048():
    game = Game2048()
    
    while not game.game_over and not game.won:
        game.display_board()
        print("Enter move (W/A/S/D) or 'q' to quit: ", end="")
        
        move_input = input().lower()
        
        if move_input == 'q':
            print("Quitting game.")
            break
        
        direction_map = {
            'w': 'up',
            'a': 'left',
            's': 'down',
            'd': 'right'
        }
        
        direction = direction_map.get(move_input)
        
        if direction:
            # The move method handles adding new tiles and checking game over for valid moves
            game.move(direction)
        else:
            print("Invalid input. Please use W, A, S, D.")
            # No game state change for invalid input, so no new tile or game over check needed here
    
    game.display_board() # Display final board state
    if game.won:
        print("Congratulations! You reached 2048!")
    elif game.game_over:
        print("Game Over! No more moves possible.")
    print(f"Final Score: {game.score}")