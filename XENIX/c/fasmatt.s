	.text
	.globl	fsincos
	.type	fsincos, @function
fsincos:
.LFB2:
	fldl	4(%esp)
	fsincos
	movl	12(%esp), %eax
	fstpl	(%eax)
	movl	16(%esp), %eax
	fstpl	(%eax)
	ret
.LFE2:
	.size	fsincos, .-fsincos
