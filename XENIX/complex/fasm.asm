TITLE fasm

.386
.387

DGROUP GROUP CONST, CONST2, _DATA, DATA, _BSS, STACK

CONST SEGMENT WORD USE32 PUBLIC 'DATA'
CONST ENDS

CONST2 SEGMENT WORD USE32 PUBLIC 'DATA'
CONST2 ENDS

_DATA SEGMENT DWORD USE32 PUBLIC 'DATA'
_DATA ENDS

DATA SEGMENT WORD USE32 PUBLIC 'DATA'
DATA ENDS

_BSS SEGMENT WORD USE32 PUBLIC 'BSS'
_BSS ENDS

STACK SEGMENT PARA USE32 STACK 'STACK'
STACK ENDS

_TEXT SEGMENT DWORD USE32 PUBLIC 'CODE'
ASSUME CS: _TEXT, DS: DGROUP, SS: DGROUP, ES: DGROUP

fsincos_ PROC NEAR
	fld qword ptr 0x4[esp]
	fsincos
	mov eax, dword ptr 0xc[esp]
	fstp qword ptr [eax]
	mov eax, dword ptr 0x10[esp]
	fstp qword ptr [eax]
	ret 0x10
fsincos_ ENDP

_TEXT ENDS

PUBLIC fsincos_

END
