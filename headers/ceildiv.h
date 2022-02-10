#pragma once

#define cdiv(x, y) ((x + (y - 1)) / (y))
#define roundup(x, y) (cdiv((x), (y)) * (y))
